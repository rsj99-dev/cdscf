use cdscf::output::{write_analysis, OutputFormat};
use cdscf::{analyze_path, AnalyzeOptions, InputFormat, Result};
use clap::{Parser, ValueEnum};
use std::fs::File;
use std::io::{self, BufWriter};
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Parser)]
#[command(
    name = "cdscf",
    version,
    about = "只统计 GenBank 或 SnapGene 文件中 CDS 编码区的密码子频率",
    after_help = "示例:\n  cdscf genomic.gbff\n  cdscf genomic.dna --format csv --output codons.csv\n  cdscf genomic.gbff --format json --exclude-stop"
)]
struct Cli {
    /// NCBI .gbff/.gb/.gbk 或 SnapGene .dna 文件
    input: PathBuf,

    /// 输出格式
    #[arg(short, long, value_enum, default_value_t = CliOutputFormat::Table)]
    format: CliOutputFormat,

    /// 输出文件；省略时写入标准输出
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// 输入格式，默认根据文件头识别
    #[arg(long, value_enum, default_value_t = CliInputFormat::Auto)]
    input_format: CliInputFormat,

    /// 纳入标记为 /pseudo 或 /pseudogene 的 CDS
    #[arg(long)]
    include_pseudo: bool,

    /// 排除 TAA、TAG、TGA，不将其计入总数和频率分母
    #[arg(long)]
    exclude_stop: bool,
}

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
enum CliOutputFormat {
    #[default]
    Table,
    Csv,
    Json,
}

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
enum CliInputFormat {
    #[default]
    Auto,
    Genbank,
    Snapgene,
}

fn run(cli: Cli) -> Result<()> {
    let input_format = match cli.input_format {
        CliInputFormat::Auto => InputFormat::Auto,
        CliInputFormat::Genbank => InputFormat::GenBank,
        CliInputFormat::Snapgene => InputFormat::SnapGene,
    };
    let options = AnalyzeOptions {
        include_pseudo: cli.include_pseudo,
        exclude_stop_codons: cli.exclude_stop,
    };
    let analysis = analyze_path(&cli.input, input_format, options)?;
    let output_format = match cli.format {
        CliOutputFormat::Table => OutputFormat::Table,
        CliOutputFormat::Csv => OutputFormat::Csv,
        CliOutputFormat::Json => OutputFormat::Json,
    };

    match cli.output {
        Some(path) => {
            let file = File::create(&path).map_err(|error| {
                cdscf::CdscfError::new(format!("无法创建输出文件“{}”：{error}", path.display()))
            })?;
            write_analysis(BufWriter::new(file), &analysis, output_format)
        }
        None => write_analysis(
            BufWriter::new(io::stdout().lock()),
            &analysis,
            output_format,
        ),
    }
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("cdscf: {error}");
            ExitCode::FAILURE
        }
    }
}
