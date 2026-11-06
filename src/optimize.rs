use crate::error::{CdscfError, Result};
use crate::stats::{amino_acid_for_codon, Analysis};
use std::collections::{BTreeMap, BTreeSet};

/// 目的序列类型；Auto 优先识别仅含 DNA 字符的序列。
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum SequenceKind {
    #[default]
    Auto,
    Dna,
    Protein,
}

/// 一个氨基酸在宿主 CDS 中计数最高的密码子。
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CodonChoice {
    pub amino_acid: char,
    pub codon: String,
    pub count: u64,
}

/// 使用宿主最高频密码子替换后的 DNA 序列。
#[derive(Clone, Debug)]
pub struct OptimizedSequence {
    pub dna: String,
    pub kind: SequenceKind,
    pub input_residues: usize,
    pub used_codons: Vec<CodonChoice>,
}

impl OptimizedSequence {
    pub fn wrapped_sequence(&self, width: usize) -> String {
        wrap_sequence(&self.dna, width.max(1))
    }

    pub fn to_fasta(&self) -> String {
        let mut fasta = String::from(">cdscf_optimized host_max_codon\n");
        fasta.push_str(&self.wrapped_sequence(60));
        fasta.push('\n');
        fasta
    }
}

/// 按宿主 CDS 中的实际计数选择最高频密码子；同计数时按密码子字母序取前者。
pub fn optimal_codons(analysis: &Analysis) -> BTreeMap<char, CodonChoice> {
    let mut choices: BTreeMap<char, CodonChoice> = BTreeMap::new();
    for entry in analysis.stats.entries() {
        if entry.count == 0 {
            continue;
        }
        let replace = match choices.get(&entry.amino_acid) {
            None => true,
            Some(current) => {
                entry.count > current.count
                    || (entry.count == current.count && entry.codon < current.codon)
            }
        };
        if replace {
            choices.insert(
                entry.amino_acid,
                CodonChoice {
                    amino_acid: entry.amino_acid,
                    codon: entry.codon,
                    count: entry.count,
                },
            );
        }
    }
    choices
}

/// 将目的 DNA 或蛋白质序列中的每个氨基酸替换为宿主最高频密码子。
pub fn optimize_sequence(
    analysis: &Analysis,
    target: &str,
    requested_kind: SequenceKind,
) -> Result<OptimizedSequence> {
    let normalized = normalize_target(target)?;
    let kind = match requested_kind {
        SequenceKind::Auto => detect_kind(&normalized),
        kind => kind,
    };
    let choices = optimal_codons(analysis);
    let mut used_amino_acids = BTreeSet::new();
    let mut dna = String::new();

    match kind {
        SequenceKind::Dna => {
            if normalized.len() % 3 != 0 {
                return Err(CdscfError::new(format!(
                    "目的 DNA 长度为 {} bp，不是 3 的整数倍，无法按密码子读取",
                    normalized.len()
                )));
            }
            for (index, codon) in normalized.chunks_exact(3).enumerate() {
                let codon = std::str::from_utf8(codon)
                    .map_err(|_| CdscfError::new("目的 DNA 含非 ASCII 字符"))?;
                let amino_acid = amino_acid_for_codon(codon).ok_or_else(|| {
                    CdscfError::new(format!(
                        "第 {} 个密码子“{codon}”含简并或无效碱基，无法保证翻译不变",
                        index + 1
                    ))
                })?;
                let choice = choices
                    .get(&amino_acid)
                    .ok_or_else(|| missing_amino_acid_error(amino_acid, "宿主 CDS"))?;
                dna.push_str(&choice.codon);
                used_amino_acids.insert(amino_acid);
            }
        }
        SequenceKind::Protein => {
            for (index, residue) in normalized.iter().copied().enumerate() {
                let amino_acid = char::from(residue);
                if !(amino_acid.is_ascii_uppercase() || amino_acid == '*') {
                    return Err(CdscfError::new(format!(
                        "蛋白质序列第 {} 位字符“{amino_acid}”无效",
                        index + 1
                    )));
                }
                if !AMINO_ACID_LETTERS.contains(&residue) {
                    return Err(CdscfError::new(format!(
                        "蛋白质序列第 {} 位“{amino_acid}”不是 20 种标准氨基酸或终止符号 *",
                        index + 1
                    )));
                }
                let choice = choices
                    .get(&amino_acid)
                    .ok_or_else(|| missing_amino_acid_error(amino_acid, "宿主 CDS"))?;
                dna.push_str(&choice.codon);
                used_amino_acids.insert(amino_acid);
            }
        }
        SequenceKind::Auto => unreachable!("目的序列类型已在上方完成识别"),
    }

    let used_codons = choices
        .into_values()
        .filter(|choice| used_amino_acids.contains(&choice.amino_acid))
        .collect();
    Ok(OptimizedSequence {
        dna,
        kind,
        input_residues: match kind {
            SequenceKind::Dna => normalized.len() / 3,
            _ => normalized.len(),
        },
        used_codons,
    })
}

fn normalize_target(target: &str) -> Result<Vec<u8>> {
    let mut normalized = Vec::new();
    for line in target.lines() {
        let line = line.trim();
        if line.starts_with('>') || line.starts_with(';') {
            continue;
        }
        for character in line.chars() {
            if character.is_ascii_alphabetic() || character == '*' {
                normalized.push(character.to_ascii_uppercase() as u8);
            } else if character.is_whitespace() || character.is_ascii_digit() {
                continue;
            } else {
                return Err(CdscfError::new(format!(
                    "目的序列含无效字符“{character}”；请输入 FASTA 或纯序列"
                )));
            }
        }
    }
    if normalized.is_empty() {
        return Err(CdscfError::new("请输入目的 DNA 或氨基酸序列"));
    }
    Ok(normalized)
}

fn detect_kind(sequence: &[u8]) -> SequenceKind {
    if sequence
        .iter()
        .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'U' | b'N'))
    {
        SequenceKind::Dna
    } else {
        SequenceKind::Protein
    }
}

fn missing_amino_acid_error(amino_acid: char, source: &str) -> CdscfError {
    CdscfError::new(format!(
        "{source}中没有氨基酸“{amino_acid}”对应的密码子计数，无法确定宿主最高频密码子"
    ))
}

fn wrap_sequence(sequence: &str, width: usize) -> String {
    sequence
        .as_bytes()
        .chunks(width)
        .map(|chunk| String::from_utf8_lossy(chunk).into_owned())
        .collect::<Vec<_>>()
        .join("\n")
}

const AMINO_ACID_LETTERS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY*";

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{CdsFeature, Record};
    use crate::stats::{analyze_records, AnalyzeOptions};

    fn host() -> Analysis {
        let sequences: [&[u8]; 2] = [b"ATGAAAAAAAAATAA", b"AAGAAGTGA"];
        let features = sequences
            .into_iter()
            .map(|sequence| CdsFeature {
                name: "host".to_owned(),
                location: "test".to_owned(),
                sequence: sequence.to_vec(),
                codon_start: 1,
                is_pseudo: false,
            })
            .collect();
        let record = Record {
            id: "host".to_owned(),
            sequence: Vec::new(),
            cds_features: features,
        };
        analyze_records(&[record], AnalyzeOptions::default()).unwrap()
    }

    #[test]
    fn chooses_highest_count_and_breaks_ties_alphabetically() {
        let choices = optimal_codons(&host());
        assert_eq!(choices[&'K'].codon, "AAA");
        assert_eq!(choices[&'K'].count, 3);
        assert_eq!(choices[&'*'].codon, "TAA");
    }

    #[test]
    fn back_translates_protein_with_host_codons() {
        let result = optimize_sequence(&host(), "MK*", SequenceKind::Protein).unwrap();
        assert_eq!(result.dna, "ATGAAATAA");
        assert_eq!(result.input_residues, 3);
        assert_eq!(result.used_codons.len(), 3);
    }

    #[test]
    fn preserves_dna_translation_while_replacing_codons() {
        let result = optimize_sequence(&host(), "ATGAAGTGA", SequenceKind::Dna).unwrap();
        assert_eq!(result.dna, "ATGAAATAA");
    }

    #[test]
    fn reads_fasta_and_wraps_output() {
        let result =
            optimize_sequence(&host(), ">target\natgaagtga\n", SequenceKind::Auto).unwrap();
        assert_eq!(result.kind, SequenceKind::Dna);
        assert_eq!(result.wrapped_sequence(3), "ATG\nAAA\nTAA");
        assert!(result.to_fasta().starts_with(">cdscf_optimized"));
    }

    #[test]
    fn rejects_incomplete_or_ambiguous_dna() {
        assert!(optimize_sequence(&host(), "ATGA", SequenceKind::Dna).is_err());
        assert!(optimize_sequence(&host(), "ATGNNTAA", SequenceKind::Dna).is_err());
    }

    #[test]
    fn reports_amino_acids_missing_from_host() {
        let error = optimize_sequence(&host(), "W", SequenceKind::Protein).unwrap_err();
        assert!(error.to_string().contains("没有氨基酸"));
    }
}
