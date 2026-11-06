mod error;
mod genbank;
mod location;
mod model;
mod optimize;
pub mod output;
mod snapgene;
mod stats;

pub use error::{CdscfError, Result};
pub use model::{CdsFeature, Record};
pub use optimize::{
    optimal_codons, optimize_sequence, CodonChoice, OptimizedSequence, SequenceKind,
};
pub use stats::{
    amino_acid_for_codon, analyze_records, Analysis, AnalyzeOptions, CodonFrequency, CodonStats,
};

use std::fs;
use std::path::Path;

/// 输入文件格式；Auto 会识别 SnapGene 二进制头，否则按 GenBank 文本解析。
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum InputFormat {
    #[default]
    Auto,
    GenBank,
    SnapGene,
}

pub fn parse_genbank(bytes: &[u8]) -> Result<Vec<Record>> {
    genbank::parse(bytes)
}

pub fn parse_snapgene(bytes: &[u8]) -> Result<Vec<Record>> {
    snapgene::parse(bytes)
}

pub fn read_records(path: impl AsRef<Path>, format: InputFormat) -> Result<Vec<Record>> {
    let path = path.as_ref();
    let bytes = fs::read(path)
        .map_err(|error| CdscfError::new(format!("无法读取文件“{}”：{error}", path.display())))?;
    match format {
        InputFormat::GenBank => parse_genbank(&bytes),
        InputFormat::SnapGene => parse_snapgene(&bytes),
        InputFormat::Auto if snapgene::looks_like(&bytes) => parse_snapgene(&bytes),
        InputFormat::Auto => parse_genbank(&bytes),
    }
}

pub fn analyze_path(
    path: impl AsRef<Path>,
    format: InputFormat,
    options: AnalyzeOptions,
) -> Result<Analysis> {
    let records = read_records(path, format)?;
    analyze_records(&records, options)
}
