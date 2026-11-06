use super::{
    amino_acid_name, CdscfApp, AMINO_ACID_ORDER, BLUE, LIGHT_BLUE, LIGHT_GRAY, NAVY, TEXT_GRAY,
};
use cdscf::{
    analyze_path, optimize_sequence, Analysis, AnalyzeOptions, InputFormat, OptimizedSequence,
    SequenceKind,
};
use eframe::egui::{self, Color32, ComboBox, Frame, RichText, ScrollArea};
use std::path::PathBuf;
use std::sync::mpsc::{self, Receiver, TryRecvError};
use std::thread;
use std::time::Duration;

type WorkerResult = Result<(Analysis, OptimizedSequence), String>;

pub(crate) struct OptimizerState {
    pub(crate) host_path: String,
    target: String,
    kind: SequenceKind,
    host_analysis: Option<Analysis>,
    result: Option<OptimizedSequence>,
    worker: Option<Receiver<WorkerResult>>,
    pub(crate) busy: bool,
}

impl Default for OptimizerState {
    fn default() -> Self {
        Self {
            host_path: String::new(),
            target: String::new(),
            kind: SequenceKind::Auto,
            host_analysis: None,
            result: None,
            worker: None,
            busy: false,
        }
    }
}

impl CdscfApp {
    pub(crate) fn poll_optimizer(&mut self) {
        let Some(receiver) = &self.optimizer.worker else {
            return;
        };
        let result = match receiver.try_recv() {
            Ok(result) => result,
            Err(TryRecvError::Empty) => return,
            Err(TryRecvError::Disconnected) => {
                self.optimizer.worker = None;
                self.optimizer.busy = false;
                self.optimizer.result = None;
                self.status = "优化线程异常终止，请重新运行".to_owned();
                self.status_is_error = true;
                return;
            }
        };
        self.optimizer.worker = None;
        self.optimizer.busy = false;
        match result {
            Ok((analysis, optimized)) => {
                self.status = format!(
                    "优化完成：{} 个输入残基，输出 {} bp DNA",
                    optimized.input_residues,
                    optimized.dna.len()
                );
                self.status_is_error = false;
                self.optimizer.host_analysis = Some(analysis);
                self.optimizer.result = Some(optimized);
            }
            Err(error) => {
                self.status = error;
                self.status_is_error = true;
                self.optimizer.result = None;
            }
        }
    }

    pub(crate) fn optimizer_panel(&mut self, context: &egui::Context) {
        egui::CentralPanel::default()
            .frame(
                Frame::default()
                    .fill(Color32::WHITE)
                    .inner_margin(egui::Margin::symmetric(24, 18)),
            )
            .show(context, |ui| {
                ScrollArea::vertical()
                    .id_salt("optimizer_scroll")
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        ui.heading(RichText::new("宿主密码子优化").color(NAVY));
                        ui.label(
                            RichText::new(
                                "先从宿主菌 CDS 确定每个氨基酸的最高频密码子，再替换目的 DNA 或蛋白质序列。",
                            )
                            .color(TEXT_GRAY),
                        );
                        ui.add_space(16.0);

                        self.host_input(ui);
                        ui.add_space(16.0);
                        self.target_input(ui);
                        ui.add_space(14.0);

                        let button = ui.add_enabled(
                            !self.optimizer.busy,
                            egui::Button::new(
                                RichText::new(if self.optimizer.busy {
                                    "正在分析宿主并优化……"
                                } else {
                                    "分析宿主并优化序列"
                                })
                                .strong()
                                .color(Color32::WHITE),
                            )
                            .fill(BLUE)
                            .min_size(egui::vec2(ui.available_width(), 38.0)),
                        );
                        if button.clicked() {
                            self.start_optimization();
                        }
                        if self.optimizer.busy {
                            ui.horizontal(|ui| {
                                ui.spinner();
                                ui.label("正在读取宿主 CDS 并计算最高频密码子");
                            });
                        }

                        if let Some(host) = &self.optimizer.host_analysis {
                            ui.add_space(12.0);
                            Frame::default()
                                .fill(LIGHT_BLUE)
                                .inner_margin(egui::Margin::symmetric(12, 8))
                                .show(ui, |ui| {
                                    ui.label(
                                        RichText::new(format!(
                                            "宿主：{} 个记录，{} 个 CDS，{} 个有效密码子",
                                            host.record_count(),
                                            host.cds_features,
                                            host.stats.total_codons
                                        ))
                                        .color(NAVY),
                                    );
                                });
                        }

                        let mut copy_sequence = false;
                        let mut export_fasta = false;
                        if let Some(result) = &self.optimizer.result {
                            ui.add_space(18.0);
                            ui.separator();
                            ui.add_space(12.0);
                            ui.heading(RichText::new("优化后的 DNA 序列").color(NAVY));
                            ui.label(
                                RichText::new(format!(
                                    "{}；输入 {} 个残基；输出 {} bp",
                                    result_kind_label(result.kind),
                                    result.input_residues,
                                    result.dna.len()
                                ))
                                .color(TEXT_GRAY),
                            );
                            ui.add_space(8.0);
                            let mut wrapped = result.wrapped_sequence(72);
                            ui.add(
                                egui::TextEdit::multiline(&mut wrapped)
                                    .font(egui::TextStyle::Monospace)
                                    .desired_rows(9)
                                    .desired_width(ui.available_width())
                                    .interactive(false),
                            );
                            ui.add_space(8.0);
                            ui.horizontal(|ui| {
                                if ui.button("复制 DNA 序列").clicked() {
                                    copy_sequence = true;
                                }
                                if ui.button("导出 FASTA").clicked() {
                                    export_fasta = true;
                                }
                            });
                            ui.add_space(16.0);
                            self.codon_choice_table(ui, result);
                        }
                        if copy_sequence {
                            if let Some(result) = &self.optimizer.result {
                                context.copy_text(result.dna.clone());
                                self.status = "优化后的 DNA 序列已复制到剪贴板".to_owned();
                                self.status_is_error = false;
                            }
                        }
                        if export_fasta {
                            self.export_optimized_fasta();
                        }
                    });
            });
        if self.optimizer.busy {
            context.request_repaint_after(Duration::from_millis(100));
        }
    }

    fn host_input(&mut self, ui: &mut egui::Ui) {
        Frame::default()
            .fill(LIGHT_GRAY)
            .inner_margin(egui::Margin::symmetric(14, 12))
            .show(ui, |ui| {
                ui.label(RichText::new("1. 输入宿主菌的基因组").strong().color(NAVY));
                ui.label("用于统计宿主每个氨基酸的最高频密码子，支持 GenBank 和 SnapGene 文件。");
                ui.add_space(8.0);
                ui.add(
                    egui::TextEdit::singleline(&mut self.optimizer.host_path)
                        .desired_width(ui.available_width())
                        .hint_text("宿主基因组 .gbff / .gb / .gbk / .dna"),
                );
                ui.add_space(6.0);
                if ui
                    .add_enabled(!self.optimizer.busy, egui::Button::new("选择宿主基因组…"))
                    .clicked()
                {
                    self.browse_host();
                }
            });
    }

    fn target_input(&mut self, ui: &mut egui::Ui) {
        Frame::default()
            .fill(LIGHT_GRAY)
            .inner_margin(egui::Margin::symmetric(14, 12))
            .show(ui, |ui| {
                ui.label(
                    RichText::new("2. 输入目的 DNA 或氨基酸序列")
                        .strong()
                        .color(NAVY),
                );
                ui.label("DNA 将在保持翻译不变的前提下替换密码子；蛋白质序列将反向翻译为宿主优化 DNA。支持 FASTA。");
                ui.add_space(8.0);
                ui.horizontal(|ui| {
                    ui.label("序列类型：");
                    ComboBox::from_id_salt("optimizer_sequence_kind")
                        .selected_text(kind_label(self.optimizer.kind))
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut self.optimizer.kind,
                                SequenceKind::Auto,
                                kind_label(SequenceKind::Auto),
                            );
                            ui.selectable_value(
                                &mut self.optimizer.kind,
                                SequenceKind::Dna,
                                kind_label(SequenceKind::Dna),
                            );
                            ui.selectable_value(
                                &mut self.optimizer.kind,
                                SequenceKind::Protein,
                                kind_label(SequenceKind::Protein),
                            );
                        });
                });
                ui.add_space(8.0);
                ui.add(
                    egui::TextEdit::multiline(&mut self.optimizer.target)
                        .font(egui::TextStyle::Monospace)
                        .desired_rows(8)
                        .desired_width(ui.available_width())
                        .hint_text("例如：ATGAAAGCT… 或 MKAV…；也支持 >sequence 开头的 FASTA"),
                );
            });
    }

    fn codon_choice_table(&self, ui: &mut egui::Ui, result: &OptimizedSequence) {
        ui.label(
            RichText::new("本次使用的宿主最高频密码子")
                .strong()
                .color(NAVY),
        );
        ui.add_space(6.0);
        let mut choices: Vec<_> = result.used_codons.iter().collect();
        choices.sort_by_key(|choice| {
            AMINO_ACID_ORDER
                .iter()
                .position(|amino_acid| *amino_acid == choice.amino_acid)
                .unwrap_or(usize::MAX)
        });
        egui::Grid::new("optimizer_codon_choices")
            .striped(true)
            .num_columns(4)
            .spacing([24.0, 6.0])
            .show(ui, |ui| {
                for heading in ["氨基酸", "名称", "宿主最高频密码子", "宿主 CDS 计数"]
                {
                    ui.label(RichText::new(heading).strong().color(NAVY));
                }
                ui.end_row();
                for choice in choices {
                    ui.label(choice.amino_acid.to_string());
                    ui.label(amino_acid_name(choice.amino_acid));
                    ui.label(RichText::new(&choice.codon).monospace().strong());
                    ui.label(choice.count.to_string());
                    ui.end_row();
                }
            });
    }

    fn browse_host(&mut self) {
        let file = rfd::FileDialog::new()
            .set_title("选择宿主菌基因组")
            .add_filter("支持的文件", &["gbff", "gb", "gbk", "dna"])
            .add_filter("GenBank", &["gbff", "gb", "gbk"])
            .add_filter("SnapGene DNA", &["dna"])
            .pick_file();
        if let Some(path) = file {
            self.optimizer.host_path = path.to_string_lossy().into_owned();
            self.status = "宿主基因组已载入，请输入目的序列并开始优化".to_owned();
            self.status_is_error = false;
        }
    }

    fn start_optimization(&mut self) {
        let path = self.optimizer.host_path.trim();
        if path.is_empty() {
            self.status = "请先选择宿主菌基因组".to_owned();
            self.status_is_error = true;
            return;
        }
        if self.optimizer.target.trim().is_empty() {
            self.status = "请输入目的 DNA 或氨基酸序列".to_owned();
            self.status_is_error = true;
            return;
        }
        let path = PathBuf::from(path);
        if !path.is_file() {
            self.status = format!("宿主基因组文件不存在：{}", path.display());
            self.status_is_error = true;
            return;
        }
        let target = self.optimizer.target.clone();
        let kind = self.optimizer.kind;
        let (sender, receiver) = mpsc::channel();
        thread::spawn(move || {
            let result = analyze_path(path, InputFormat::Auto, AnalyzeOptions::default())
                .map_err(|error| error.to_string())
                .and_then(|analysis| {
                    optimize_sequence(&analysis, &target, kind)
                        .map(|optimized| (analysis, optimized))
                        .map_err(|error| error.to_string())
                });
            let _ = sender.send(result);
        });
        self.optimizer.worker = Some(receiver);
        self.optimizer.busy = true;
        self.optimizer.result = None;
        self.status = "正在分析宿主 CDS 并优化目的序列……".to_owned();
        self.status_is_error = false;
    }

    fn export_optimized_fasta(&mut self) {
        let Some(result) = &self.optimizer.result else {
            return;
        };
        let destination = rfd::FileDialog::new()
            .set_title("导出优化后的 DNA")
            .set_file_name("cdscf_optimized.fasta")
            .add_filter("FASTA 序列", &["fasta", "fa"])
            .save_file();
        let Some(destination) = destination else {
            return;
        };
        match std::fs::write(&destination, result.to_fasta()) {
            Ok(()) => {
                self.status = format!("优化序列已导出至 {}", destination.display());
                self.status_is_error = false;
            }
            Err(error) => {
                self.status = format!("无法导出 FASTA：{error}");
                self.status_is_error = true;
            }
        }
    }
}

fn kind_label(kind: SequenceKind) -> &'static str {
    match kind {
        SequenceKind::Auto => "自动识别",
        SequenceKind::Dna => "DNA",
        SequenceKind::Protein => "蛋白质 / 氨基酸",
    }
}

fn result_kind_label(kind: SequenceKind) -> &'static str {
    match kind {
        SequenceKind::Dna => "DNA 密码子替换",
        SequenceKind::Protein => "蛋白质反向翻译",
        SequenceKind::Auto => "序列优化",
    }
}
