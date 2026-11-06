#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod optimizer_gui;

use cdscf::output::{write_analysis, OutputFormat};
use cdscf::{analyze_path, Analysis, AnalyzeOptions, CodonFrequency, InputFormat};
use eframe::egui::{
    self, Color32, ComboBox, FontData, FontDefinitions, FontFamily, Frame, ProgressBar, RichText,
    ScrollArea, Stroke, TextureHandle, TextureOptions,
};
use optimizer_gui::OptimizerState;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::sync::mpsc::{self, Receiver, TryRecvError};
use std::thread;
use std::time::Duration;

const NAVY: Color32 = Color32::from_rgb(15, 48, 87);
const BLUE: Color32 = Color32::from_rgb(17, 119, 187);
const LIGHT_BLUE: Color32 = Color32::from_rgb(232, 244, 252);
const LIGHT_GRAY: Color32 = Color32::from_rgb(246, 248, 251);
const TEXT_GRAY: Color32 = Color32::from_rgb(83, 96, 110);
const ERROR_RED: Color32 = Color32::from_rgb(180, 45, 45);
const AMINO_ACID_ORDER: [char; 21] = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y',
    'V', '*',
];

fn main() -> eframe::Result<()> {
    let icon = eframe::icon_data::from_png_bytes(include_bytes!("../assets/cdscf_icon.png"))
        .expect("内置应用图标应可解码");
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("cdscf — CDS Codon Frequency")
            .with_inner_size([1180.0, 780.0])
            .with_min_inner_size([960.0, 650.0])
            .with_icon(icon),
        ..Default::default()
    };
    eframe::run_native(
        "cdscf",
        options,
        Box::new(|context| Ok(Box::new(CdscfApp::new(context)))),
    )
}

#[derive(Clone, Copy, Default, Eq, PartialEq)]
enum InputChoice {
    #[default]
    Auto,
    GenBank,
    SnapGene,
}

impl InputChoice {
    fn label(self) -> &'static str {
        match self {
            Self::Auto => "自动识别",
            Self::GenBank => "GenBank / GBFF",
            Self::SnapGene => "SnapGene DNA",
        }
    }

    fn format(self) -> InputFormat {
        match self {
            Self::Auto => InputFormat::Auto,
            Self::GenBank => InputFormat::GenBank,
            Self::SnapGene => InputFormat::SnapGene,
        }
    }
}

#[derive(Clone, Copy, Default, Eq, PartialEq)]
enum OutputChoice {
    #[default]
    Csv,
    Json,
    Table,
}

impl OutputChoice {
    fn label(self) -> &'static str {
        match self {
            Self::Csv => "CSV 表格",
            Self::Json => "JSON 数据",
            Self::Table => "文本表格",
        }
    }

    fn extension(self) -> &'static str {
        match self {
            Self::Csv => "csv",
            Self::Json => "json",
            Self::Table => "txt",
        }
    }

    fn format(self) -> OutputFormat {
        match self {
            Self::Csv => OutputFormat::Csv,
            Self::Json => OutputFormat::Json,
            Self::Table => OutputFormat::Table,
        }
    }
}

#[derive(Clone, Copy, Default, Eq, PartialEq)]
enum SortChoice {
    #[default]
    AminoAcidFrequency,
    Codon,
    CountDescending,
    FrequencyDescending,
}

impl SortChoice {
    fn label(self) -> &'static str {
        match self {
            Self::AminoAcidFrequency => "按氨基酸分组（高频优先）",
            Self::Codon => "按密码子排序",
            Self::CountDescending => "按计数降序",
            Self::FrequencyDescending => "按频率降序",
        }
    }
}

#[derive(Clone, Copy, Default, Eq, PartialEq)]
enum WorkspaceTab {
    #[default]
    Frequency,
    Optimizer,
}

struct CdscfApp {
    logo: TextureHandle,
    active_tab: WorkspaceTab,
    optimizer: OptimizerState,
    input_path: String,
    input_choice: InputChoice,
    output_choice: OutputChoice,
    include_pseudo: bool,
    exclude_stop: bool,
    analysis: Option<Analysis>,
    worker: Option<Receiver<Result<Analysis, String>>>,
    busy: bool,
    status: String,
    status_is_error: bool,
    filter: String,
    sort: SortChoice,
    show_about: bool,
}

impl CdscfApp {
    fn new(context: &eframe::CreationContext<'_>) -> Self {
        context.egui_ctx.set_visuals(egui::Visuals::light());
        install_chinese_font(&context.egui_ctx);
        let logo = load_logo(&context.egui_ctx);
        Self {
            logo,
            active_tab: WorkspaceTab::Frequency,
            optimizer: OptimizerState::default(),
            input_path: String::new(),
            input_choice: InputChoice::Auto,
            output_choice: OutputChoice::Csv,
            include_pseudo: false,
            exclude_stop: false,
            analysis: None,
            worker: None,
            busy: false,
            status: "请选择 GenBank 或 SnapGene 文件".to_owned(),
            status_is_error: false,
            filter: String::new(),
            sort: SortChoice::AminoAcidFrequency,
            show_about: false,
        }
    }

    fn poll_worker(&mut self) {
        let Some(receiver) = &self.worker else {
            return;
        };
        let result = match receiver.try_recv() {
            Ok(result) => result,
            Err(TryRecvError::Empty) => return,
            Err(TryRecvError::Disconnected) => {
                self.worker = None;
                self.busy = false;
                self.analysis = None;
                self.status = "分析线程异常终止，请重新运行分析".to_owned();
                self.status_is_error = true;
                return;
            }
        };
        self.worker = None;
        self.busy = false;
        match result {
            Ok(analysis) => {
                self.status = format!(
                    "分析完成：{} 个 CDS，{} 个有效密码子",
                    analysis.cds_features, analysis.stats.total_codons
                );
                self.status_is_error = false;
                self.analysis = Some(analysis);
            }
            Err(error) => {
                self.status = error;
                self.status_is_error = true;
                self.analysis = None;
            }
        }
    }

    fn handle_dropped_file(&mut self, context: &egui::Context) {
        let dropped = context.input(|input| {
            input
                .raw
                .dropped_files
                .first()
                .and_then(|file| file.path.clone())
        });
        let Some(path) = dropped else {
            return;
        };
        match self.active_tab {
            WorkspaceTab::Frequency if !self.busy => {
                self.input_path = path.to_string_lossy().into_owned();
                self.status = "文件已载入，点击“开始分析”".to_owned();
            }
            WorkspaceTab::Optimizer if !self.optimizer.busy => {
                self.optimizer.host_path = path.to_string_lossy().into_owned();
                self.status = "宿主基因组已载入，请输入目的序列并开始优化".to_owned();
            }
            _ => return,
        }
        self.status_is_error = false;
    }

    fn browse_input(&mut self) {
        let file = rfd::FileDialog::new()
            .set_title("选择基因组或 DNA 文件")
            .add_filter("支持的文件", &["gbff", "gb", "gbk", "dna"])
            .add_filter("GenBank", &["gbff", "gb", "gbk"])
            .add_filter("SnapGene DNA", &["dna"])
            .pick_file();
        if let Some(path) = file {
            self.input_path = path.to_string_lossy().into_owned();
            self.status = "文件已载入，点击“开始分析”".to_owned();
            self.status_is_error = false;
        }
    }

    fn start_analysis(&mut self) {
        let path = self.input_path.trim();
        if path.is_empty() {
            self.status = "请先选择输入文件".to_owned();
            self.status_is_error = true;
            return;
        }
        let path = PathBuf::from(path);
        if !path.is_file() {
            self.status = format!("文件不存在：{}", path.display());
            self.status_is_error = true;
            return;
        }

        let (sender, receiver) = mpsc::channel();
        let input_format = self.input_choice.format();
        let options = AnalyzeOptions {
            include_pseudo: self.include_pseudo,
            exclude_stop_codons: self.exclude_stop,
        };
        thread::spawn(move || {
            let result =
                analyze_path(path, input_format, options).map_err(|error| error.to_string());
            let _ = sender.send(result);
        });
        self.worker = Some(receiver);
        self.busy = true;
        self.analysis = None;
        self.status = "正在读取文件并统计 CDS 密码子……".to_owned();
        self.status_is_error = false;
    }

    fn export_result(&mut self) {
        let Some(analysis) = &self.analysis else {
            self.status = "没有可导出的分析结果".to_owned();
            self.status_is_error = true;
            return;
        };
        let extension = self.output_choice.extension();
        let file_name = format!("codon_frequency.{extension}");
        let destination = rfd::FileDialog::new()
            .set_title("导出密码子频率")
            .set_file_name(&file_name)
            .add_filter(self.output_choice.label(), &[extension])
            .save_file();
        let Some(destination) = destination else {
            return;
        };

        let result = File::create(&destination)
            .map_err(|error| format!("无法创建输出文件：{error}"))
            .and_then(|file| {
                write_analysis(BufWriter::new(file), analysis, self.output_choice.format())
                    .map_err(|error| error.to_string())
            });
        match result {
            Ok(()) => {
                self.status = format!("结果已导出至 {}", destination.display());
                self.status_is_error = false;
            }
            Err(error) => {
                self.status = error;
                self.status_is_error = true;
            }
        }
    }

    fn header(&mut self, context: &egui::Context) {
        egui::TopBottomPanel::top("header")
            .exact_height(124.0)
            .frame(
                Frame::default()
                    .fill(Color32::WHITE)
                    .inner_margin(egui::Margin::symmetric(20, 8))
                    .stroke(Stroke::new(1.0_f32, Color32::from_rgb(225, 231, 237))),
            )
            .show(context, |ui| {
                ui.horizontal(|ui| {
                    ui.add(
                        egui::Image::new(&self.logo).fit_to_exact_size(egui::vec2(190.0, 106.0)),
                    );
                    ui.add_space(20.0);
                    ui.vertical(|ui| {
                        ui.add_space(20.0);
                        ui.label(
                            RichText::new("CDS Codon Frequency")
                                .size(28.0)
                                .strong()
                                .color(NAVY),
                        );
                        ui.label(
                            RichText::new("只从 GenBank 与 SnapGene 注释的编码区统计密码子频率")
                                .size(14.0)
                                .color(TEXT_GRAY),
                        );
                    });
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.button("关于").clicked() {
                            self.show_about = true;
                        }
                    });
                });
            });
    }

    fn tab_bar(&mut self, context: &egui::Context) {
        egui::TopBottomPanel::top("workspace_tabs")
            .exact_height(42.0)
            .frame(Frame::default().fill(Color32::WHITE))
            .show(context, |ui| {
                ui.horizontal(|ui| {
                    ui.add_space(18.0);
                    ui.selectable_value(
                        &mut self.active_tab,
                        WorkspaceTab::Frequency,
                        RichText::new("CDS 密码子频率").size(16.0).color(NAVY),
                    );
                    ui.add_space(12.0);
                    ui.selectable_value(
                        &mut self.active_tab,
                        WorkspaceTab::Optimizer,
                        RichText::new("宿主密码子优化").size(16.0).color(NAVY),
                    );
                });
            });
    }

    fn controls(&mut self, context: &egui::Context) {
        egui::SidePanel::left("controls")
            .exact_width(310.0)
            .frame(
                Frame::default()
                    .fill(LIGHT_GRAY)
                    .inner_margin(egui::Margin::symmetric(18, 18)),
            )
            .show(context, |ui| {
                ui.heading(RichText::new("分析设置").color(NAVY));
                ui.add_space(16.0);

                ui.label(RichText::new("输入文件").strong());
                ui.add(
                    egui::TextEdit::singleline(&mut self.input_path)
                        .desired_width(274.0)
                        .hint_text("选择或拖入 .gbff / .dna 文件"),
                );
                ui.add_space(6.0);
                if ui
                    .add_enabled(!self.busy, egui::Button::new("浏览文件…").min_size(egui::vec2(274.0, 30.0)))
                    .clicked()
                {
                    self.browse_input();
                }
                ui.label(RichText::new("也可以将文件直接拖放到窗口中").small().color(TEXT_GRAY));
                ui.add_space(14.0);

                ui.label(RichText::new("输入格式").strong());
                ComboBox::from_id_salt("input_format")
                    .selected_text(self.input_choice.label())
                    .width(274.0)
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut self.input_choice, InputChoice::Auto, InputChoice::Auto.label());
                        ui.selectable_value(&mut self.input_choice, InputChoice::GenBank, InputChoice::GenBank.label());
                        ui.selectable_value(&mut self.input_choice, InputChoice::SnapGene, InputChoice::SnapGene.label());
                    });
                ui.add_space(14.0);

                ui.checkbox(&mut self.include_pseudo, "纳入伪基因 CDS");
                ui.checkbox(&mut self.exclude_stop, "排除终止密码子 TAA / TAG / TGA");
                ui.add_space(18.0);

                let analyze = ui.add_enabled(
                    !self.busy,
                    egui::Button::new(
                        RichText::new(if self.busy { "正在分析……" } else { "开始分析" })
                            .size(17.0)
                            .strong()
                            .color(Color32::WHITE),
                    )
                    .fill(BLUE)
                    .min_size(egui::vec2(274.0, 42.0)),
                );
                if analyze.clicked() {
                    self.start_analysis();
                }
                if self.busy {
                    ui.add_space(10.0);
                    ui.horizontal(|ui| {
                        ui.spinner();
                        ui.label("大文件分析在后台运行，界面仍可操作");
                    });
                }

                ui.add_space(22.0);
                ui.separator();
                ui.add_space(14.0);
                ui.label(RichText::new("导出格式").strong());
                ComboBox::from_id_salt("output_format")
                    .selected_text(self.output_choice.label())
                    .width(274.0)
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut self.output_choice, OutputChoice::Csv, OutputChoice::Csv.label());
                        ui.selectable_value(&mut self.output_choice, OutputChoice::Json, OutputChoice::Json.label());
                        ui.selectable_value(&mut self.output_choice, OutputChoice::Table, OutputChoice::Table.label());
                    });
                ui.add_space(8.0);
                if ui
                    .add_enabled(
                        self.analysis.is_some() && !self.busy,
                        egui::Button::new("导出分析结果").min_size(egui::vec2(274.0, 34.0)),
                    )
                    .clicked()
                {
                    self.export_result();
                }

                ui.add_space(20.0);
                ui.separator();
                ui.add_space(10.0);
                ui.label(RichText::new("统计口径").strong().color(NAVY));
                ui.label(
                    RichText::new(
                        "仅统计完整 ACGT 三联体；模糊密码子和末端不足三个碱基的部分不计入频率分母。氨基酸分类采用标准遗传密码表。",
                    )
                    .small()
                    .color(TEXT_GRAY),
                );
            });
    }

    fn results(&mut self, context: &egui::Context) {
        egui::CentralPanel::default()
            .frame(
                            Frame::default()
                                .fill(Color32::WHITE)
                                .inner_margin(egui::Margin::symmetric(22, 18)),
                        )
            .show(context, |ui| {
                let Some(analysis) = self.analysis.clone() else {
                    self.empty_results(ui);
                    return;
                };

                ui.horizontal(|ui| {
                    ui.heading(RichText::new("密码子频率").color(NAVY));
                    ui.label(RichText::new(format!("记录：{}", analysis.record_ids.join(", "))).color(TEXT_GRAY));
                });
                ui.add_space(12.0);
                summary_cards(ui, &analysis);
                ui.add_space(12.0);
                ui.label(
                    RichText::new(format!(
                        "跳过模糊密码子 {} 个；末端余碱基 {} 个；排除终止密码子 {} 个；跳过伪基因 CDS {} 个",
                        analysis.stats.ambiguous_codons,
                        analysis.stats.trailing_bases,
                        analysis.stats.excluded_stop_codons,
                        analysis.skipped_pseudo
                    ))
                    .small()
                    .color(TEXT_GRAY),
                );
                ui.add_space(14.0);

                ui.horizontal(|ui| {
                    ui.label("筛选密码子：");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.filter)
                            .desired_width(110.0)
                            .hint_text("例如 ATG"),
                    );
                    ui.add_space(12.0);
                    ComboBox::from_id_salt("sort_order")
                        .selected_text(self.sort.label())
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut self.sort, SortChoice::AminoAcidFrequency, SortChoice::AminoAcidFrequency.label());
                            ui.selectable_value(&mut self.sort, SortChoice::Codon, SortChoice::Codon.label());
                            ui.selectable_value(&mut self.sort, SortChoice::CountDescending, SortChoice::CountDescending.label());
                            ui.selectable_value(&mut self.sort, SortChoice::FrequencyDescending, SortChoice::FrequencyDescending.label());
                        });
                });
                ui.add_space(8.0);
                self.codon_table(ui, &analysis);
            });
    }

    fn empty_results(&self, ui: &mut egui::Ui) {
        ui.vertical_centered(|ui| {
            ui.add_space(150.0);
            ui.label(
                RichText::new("cdscf")
                    .size(42.0)
                    .strong()
                    .color(Color32::from_rgb(200, 215, 227)),
            );
            ui.add_space(10.0);
            ui.label(
                RichText::new("选择输入文件并开始分析")
                    .size(22.0)
                    .color(NAVY),
            );
            ui.label(RichText::new("结果表将在这里显示全部 64 种密码子").color(TEXT_GRAY));
        });
    }

    fn codon_table(&mut self, ui: &mut egui::Ui, analysis: &Analysis) {
        let query = self.filter.trim().to_ascii_uppercase();
        let entries: Vec<_> = analysis
            .stats
            .entries()
            .filter(|entry| query.is_empty() || entry.codon.contains(&query))
            .collect();
        let maximum = entries
            .iter()
            .map(|entry| entry.frequency)
            .fold(0.0_f64, f64::max)
            .max(f64::MIN_POSITIVE);

        ScrollArea::vertical()
            .id_salt("codon_results")
            .auto_shrink([false, false])
            .show(ui, |ui| {
                if self.sort == SortChoice::AminoAcidFrequency {
                    amino_acid_groups(ui, &entries, maximum);
                } else {
                    let mut entries = entries;
                    entries.sort_by(|left, right| match self.sort {
                        SortChoice::Codon | SortChoice::AminoAcidFrequency => {
                            left.codon.cmp(&right.codon)
                        }
                        SortChoice::CountDescending => right
                            .count
                            .cmp(&left.count)
                            .then_with(|| left.codon.cmp(&right.codon)),
                        SortChoice::FrequencyDescending => right
                            .frequency
                            .total_cmp(&left.frequency)
                            .then_with(|| left.codon.cmp(&right.codon)),
                    });
                    codon_grid(ui, "codon_grid", entries.iter(), maximum, true);
                }
            });
    }

    fn status_bar(&self, context: &egui::Context) {
        egui::TopBottomPanel::bottom("status")
            .exact_height(30.0)
            .frame(Frame::default().fill(Color32::WHITE))
            .show(context, |ui| {
                ui.horizontal_centered(|ui| {
                    let color = if self.status_is_error {
                        ERROR_RED
                    } else {
                        TEXT_GRAY
                    };
                    ui.label(RichText::new(&self.status).small().color(color));
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        ui.label(
                            RichText::new(format!("v{}", env!("CARGO_PKG_VERSION")))
                                .small()
                                .color(TEXT_GRAY),
                        );
                    });
                });
            });
    }

    fn about(&mut self, context: &egui::Context) {
        if !self.show_about {
            return;
        }
        egui::Window::new("关于 cdscf")
            .collapsible(false)
            .resizable(false)
            .default_width(390.0)
            .open(&mut self.show_about)
            .show(context, |ui| {
                ui.vertical_centered(|ui| {
                    ui.add(
                        egui::Image::new(&self.logo).fit_to_exact_size(egui::vec2(260.0, 146.0)),
                    );
                    ui.label(RichText::new("cdscf").size(26.0).strong().color(NAVY));
                    ui.label(format!("版本 {}", env!("CARGO_PKG_VERSION")));
                    ui.add_space(8.0);
                    ui.label("CDS Codon Frequency");
                    ui.label("Rust / eframe / egui");
                    ui.add_space(8.0);
                    ui.label("仅统计注释编码区中的完整标准密码子");
                    ui.label("氨基酸分类采用标准遗传密码表");
                    ui.label("许可证：MIT");
                    ui.add_space(12.0);
                    ui.separator();
                    ui.add_space(8.0);
                    ui.label(
                        RichText::new("南京工业大学生物能源研究所")
                            .strong()
                            .color(NAVY),
                    );
                    ui.label("周俊课题组");
                    ui.label("合成生物与生物信息学小组");
                    ui.add_space(6.0);
                    ui.hyperlink_to(
                        "问题汇报邮箱：rsj1999@njtech.edu.cn",
                        "mailto:rsj1999@njtech.edu.cn",
                    );
                });
            });
    }
}

fn amino_acid_groups(ui: &mut egui::Ui, entries: &[CodonFrequency], maximum: f64) {
    for amino_acid in AMINO_ACID_ORDER {
        let mut group: Vec<_> = entries
            .iter()
            .filter(|entry| entry.amino_acid == amino_acid)
            .collect();
        if group.is_empty() {
            continue;
        }
        group.sort_by(|left, right| {
            right
                .frequency
                .total_cmp(&left.frequency)
                .then_with(|| left.codon.cmp(&right.codon))
        });
        let title = format!(
            "{} · {}（{} 个密码子）",
            amino_acid,
            amino_acid_name(amino_acid),
            group.len()
        );
        egui::CollapsingHeader::new(RichText::new(title).strong().color(NAVY))
            .id_salt(format!("amino_acid_header_{amino_acid}"))
            .default_open(true)
            .show(ui, |ui| {
                codon_grid(
                    ui,
                    format!("amino_acid_grid_{amino_acid}"),
                    group,
                    maximum,
                    false,
                );
            });
        ui.add_space(6.0);
    }
}

fn codon_grid<'a>(
    ui: &mut egui::Ui,
    id: impl std::hash::Hash,
    entries: impl IntoIterator<Item = &'a CodonFrequency>,
    maximum: f64,
    show_amino_acid: bool,
) {
    egui::Grid::new(id)
        .striped(true)
        .num_columns(if show_amino_acid { 6 } else { 5 })
        .min_col_width(68.0)
        .spacing([18.0, 7.0])
        .show(ui, |ui| {
            if show_amino_acid {
                for heading in [
                    "密码子",
                    "氨基酸",
                    "计数",
                    "相对频率",
                    "每千密码子",
                    "相对强度",
                ] {
                    ui.label(RichText::new(heading).strong().color(NAVY));
                }
            } else {
                for heading in ["密码子", "计数", "相对频率", "每千密码子", "相对强度"]
                {
                    ui.label(RichText::new(heading).strong().color(NAVY));
                }
            }
            ui.end_row();
            for entry in entries {
                ui.label(RichText::new(&entry.codon).monospace().strong());
                if show_amino_acid {
                    ui.label(format!(
                        "{} {}",
                        entry.amino_acid,
                        amino_acid_name(entry.amino_acid)
                    ));
                }
                ui.label(entry.count.to_string());
                ui.label(format!("{:.8}", entry.frequency));
                ui.label(format!("{:.3}", entry.frequency * 1000.0));
                ui.add(
                    ProgressBar::new((entry.frequency / maximum) as f32)
                        .desired_width(150.0)
                        .fill(BLUE),
                );
                ui.end_row();
            }
        });
}

fn amino_acid_name(amino_acid: char) -> &'static str {
    match amino_acid {
        'A' => "丙氨酸",
        'R' => "精氨酸",
        'N' => "天冬酰胺",
        'D' => "天冬氨酸",
        'C' => "半胱氨酸",
        'Q' => "谷氨酰胺",
        'E' => "谷氨酸",
        'G' => "甘氨酸",
        'H' => "组氨酸",
        'I' => "异亮氨酸",
        'L' => "亮氨酸",
        'K' => "赖氨酸",
        'M' => "甲硫氨酸",
        'F' => "苯丙氨酸",
        'P' => "脯氨酸",
        'S' => "丝氨酸",
        'T' => "苏氨酸",
        'W' => "色氨酸",
        'Y' => "酪氨酸",
        'V' => "缬氨酸",
        '*' => "终止密码子",
        _ => "未知",
    }
}

impl eframe::App for CdscfApp {
    fn update(&mut self, context: &egui::Context, _frame: &mut eframe::Frame) {
        self.poll_worker();
        self.poll_optimizer();
        self.handle_dropped_file(context);
        self.header(context);
        self.tab_bar(context);
        self.status_bar(context);
        match self.active_tab {
            WorkspaceTab::Frequency => {
                self.controls(context);
                self.results(context);
            }
            WorkspaceTab::Optimizer => self.optimizer_panel(context),
        }
        self.about(context);
        if self.busy {
            context.request_repaint_after(Duration::from_millis(100));
        }
    }
}

fn summary_cards(ui: &mut egui::Ui, analysis: &Analysis) {
    let cards = [
        ("记录数", analysis.record_count().to_string()),
        ("CDS 数量", analysis.cds_features.to_string()),
        ("有效密码子", analysis.stats.total_codons.to_string()),
        ("模糊密码子", analysis.stats.ambiguous_codons.to_string()),
    ];
    ui.horizontal_wrapped(|ui| {
        for (label, value) in cards {
            Frame::default()
                .fill(LIGHT_BLUE)
                .inner_margin(egui::Margin::symmetric(16, 10))
                .show(ui, |ui| {
                    ui.set_min_size(egui::vec2(145.0, 52.0));
                    ui.vertical(|ui| {
                        ui.label(RichText::new(label).small().color(TEXT_GRAY));
                        ui.label(RichText::new(value).size(21.0).strong().color(NAVY));
                    });
                });
            ui.add_space(8.0);
        }
    });
}

fn load_logo(context: &egui::Context) -> TextureHandle {
    let image = image::load_from_memory(include_bytes!("../assets/cdscf_logo.png"))
        .expect("内置 Logo 应可解码")
        .to_rgba8();
    let size = [image.width() as usize, image.height() as usize];
    let image = egui::ColorImage::from_rgba_unmultiplied(size, image.as_raw());
    context.load_texture("cdscf_logo", image, TextureOptions::LINEAR)
}

fn install_chinese_font(context: &egui::Context) {
    let candidates = [
        Path::new(r"C:\Windows\Fonts\msyh.ttc"),
        Path::new(r"C:\Windows\Fonts\simhei.ttf"),
        Path::new(r"C:\Windows\Fonts\simsun.ttc"),
    ];
    let Some(bytes) = candidates.iter().find_map(|path| std::fs::read(path).ok()) else {
        return;
    };
    let mut fonts = FontDefinitions::default();
    fonts
        .font_data
        .insert("windows_cjk".to_owned(), FontData::from_owned(bytes).into());
    for family in [FontFamily::Proportional, FontFamily::Monospace] {
        fonts
            .families
            .entry(family)
            .or_default()
            .insert(0, "windows_cjk".to_owned());
    }
    context.set_fonts(fonts);
}
