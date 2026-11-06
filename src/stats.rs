use crate::error::{CdscfError, Result};
use crate::model::Record;

const BASES: [u8; 4] = *b"ACGT";
// 标准遗传密码表，顺序与 AAA…TTT 的密码子索引一致。
const AMINO_ACIDS: [char; 64] = [
    'K', 'N', 'K', 'N', 'T', 'T', 'T', 'T', 'R', 'S', 'R', 'S', 'I', 'I', 'M', 'I', 'Q', 'H', 'Q',
    'H', 'P', 'P', 'P', 'P', 'R', 'R', 'R', 'R', 'L', 'L', 'L', 'L', 'E', 'D', 'E', 'D', 'A', 'A',
    'A', 'A', 'G', 'G', 'G', 'G', 'V', 'V', 'V', 'V', '*', 'Y', '*', 'Y', 'S', 'S', 'S', 'S', '*',
    'C', 'W', 'C', 'L', 'F', 'L', 'F',
];

/// 控制 CDS 筛选和密码子计数方式。
#[derive(Clone, Copy, Debug, Default)]
pub struct AnalyzeOptions {
    pub include_pseudo: bool,
    pub exclude_stop_codons: bool,
}

/// 一个标准密码子的计数与相对频率。
#[derive(Clone, Debug, PartialEq)]
pub struct CodonFrequency {
    pub codon: String,
    pub amino_acid: char,
    pub count: u64,
    pub frequency: f64,
}

/// 64 种标准 DNA 密码子的汇总计数。
#[derive(Clone, Debug)]
pub struct CodonStats {
    counts: [u64; 64],
    pub total_codons: u64,
    pub ambiguous_codons: u64,
    pub trailing_bases: u64,
    pub excluded_stop_codons: u64,
}

impl Default for CodonStats {
    fn default() -> Self {
        Self {
            counts: [0; 64],
            total_codons: 0,
            ambiguous_codons: 0,
            trailing_bases: 0,
            excluded_stop_codons: 0,
        }
    }
}

impl CodonStats {
    pub fn count(&self, codon: &str) -> Option<u64> {
        let bytes = codon.as_bytes();
        if bytes.len() != 3 {
            return None;
        }
        let normalized = [
            bytes[0].to_ascii_uppercase(),
            bytes[1].to_ascii_uppercase(),
            bytes[2].to_ascii_uppercase(),
        ];
        codon_index(&normalized).map(|index| self.counts[index])
    }

    pub fn entries(&self) -> impl Iterator<Item = CodonFrequency> + '_ {
        self.counts
            .iter()
            .enumerate()
            .map(|(index, count)| CodonFrequency {
                codon: codon_name(index),
                amino_acid: AMINO_ACIDS[index],
                count: *count,
                frequency: if self.total_codons == 0 {
                    0.0
                } else {
                    *count as f64 / self.total_codons as f64
                },
            })
    }

    fn add_cds(&mut self, sequence: &[u8], codon_start: usize, exclude_stops: bool) {
        let offset = codon_start - 1;
        if offset >= sequence.len() {
            self.trailing_bases += sequence.len() as u64;
            return;
        }
        let coding_sequence = &sequence[offset..];
        self.trailing_bases += (coding_sequence.len() % 3) as u64;
        for codon in coding_sequence.chunks_exact(3) {
            let Some(index) = codon_index(codon) else {
                self.ambiguous_codons += 1;
                continue;
            };
            if exclude_stops
                && matches!(
                    codon,
                    [b'T', b'A', b'A'] | [b'T', b'A', b'G'] | [b'T', b'G', b'A']
                )
            {
                self.excluded_stop_codons += 1;
                continue;
            }
            self.counts[index] += 1;
            self.total_codons += 1;
        }
    }
}

/// 所有输入记录的汇总分析结果。
#[derive(Clone, Debug)]
pub struct Analysis {
    pub record_ids: Vec<String>,
    pub cds_features: usize,
    pub skipped_pseudo: usize,
    pub stats: CodonStats,
}

impl Analysis {
    pub fn record_count(&self) -> usize {
        self.record_ids.len()
    }
}

pub fn analyze_records(records: &[Record], options: AnalyzeOptions) -> Result<Analysis> {
    let mut analysis = Analysis {
        record_ids: records.iter().map(|record| record.id.clone()).collect(),
        cds_features: 0,
        skipped_pseudo: 0,
        stats: CodonStats::default(),
    };

    for record in records {
        for feature in &record.cds_features {
            if feature.is_pseudo && !options.include_pseudo {
                analysis.skipped_pseudo += 1;
                continue;
            }
            if !(1..=3).contains(&feature.codon_start) {
                return Err(CdscfError::new(format!(
                    "记录 {} 的 CDS“{}”具有无效 codon_start={}",
                    record.id, feature.name, feature.codon_start
                )));
            }
            analysis.cds_features += 1;
            analysis.stats.add_cds(
                &feature.sequence,
                feature.codon_start,
                options.exclude_stop_codons,
            );
        }
    }
    if analysis.cds_features == 0 {
        let hint = if analysis.skipped_pseudo > 0 {
            "；可使用 include_pseudo 选项纳入伪基因 CDS"
        } else {
            ""
        };
        return Err(CdscfError::new(format!(
            "输入文件中没有可统计的 CDS 编码区{hint}"
        )));
    }
    Ok(analysis)
}

/// 使用标准遗传密码表返回密码子对应的氨基酸，终止密码子返回 `*`。
pub fn amino_acid_for_codon(codon: &str) -> Option<char> {
    let bytes = codon.as_bytes();
    if bytes.len() != 3 {
        return None;
    }
    let normalized = [
        bytes[0].to_ascii_uppercase(),
        bytes[1].to_ascii_uppercase(),
        bytes[2].to_ascii_uppercase(),
    ];
    codon_index(&normalized).map(|index| AMINO_ACIDS[index])
}

fn codon_index(codon: &[u8]) -> Option<usize> {
    if codon.len() != 3 {
        return None;
    }
    let mut index = 0usize;
    for base in codon {
        let value = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        index = index * 4 + value;
    }
    Some(index)
}

fn codon_name(mut index: usize) -> String {
    let mut codon = [b'A'; 3];
    for position in (0..3).rev() {
        codon[position] = BASES[index % 4];
        index /= 4;
    }
    String::from_utf8(codon.to_vec()).unwrap_or_else(|_| "NNN".to_owned())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::CdsFeature;

    fn record(sequence: &[u8], codon_start: usize, pseudo: bool) -> Record {
        Record {
            id: "test".to_owned(),
            sequence: sequence.to_vec(),
            cds_features: vec![CdsFeature {
                name: "cds".to_owned(),
                location: "test".to_owned(),
                sequence: sequence.to_vec(),
                codon_start,
                is_pseudo: pseudo,
            }],
        }
    }

    #[test]
    fn counts_complete_standard_codons_only() {
        let records = vec![record(b"ATGAAATAANNNCC", 1, false)];
        let analysis = analyze_records(&records, AnalyzeOptions::default()).unwrap();
        assert_eq!(analysis.stats.total_codons, 3);
        assert_eq!(analysis.stats.ambiguous_codons, 1);
        assert_eq!(analysis.stats.trailing_bases, 2);
        assert_eq!(analysis.stats.count("ATG"), Some(1));
        assert_eq!(analysis.stats.count("AAA"), Some(1));
        assert_eq!(analysis.stats.count("TAA"), Some(1));
        assert_eq!(analysis.stats.count("atg"), Some(1));
    }

    #[test]
    fn honors_codon_start() {
        let records = vec![record(b"AATGAAA", 2, false)];
        let analysis = analyze_records(&records, AnalyzeOptions::default()).unwrap();
        assert_eq!(analysis.stats.count("ATG"), Some(1));
        assert_eq!(analysis.stats.count("AAA"), Some(1));
        assert_eq!(analysis.stats.trailing_bases, 0);
    }

    #[test]
    fn skips_pseudo_unless_requested() {
        let records = vec![record(b"ATG", 1, true)];
        assert!(analyze_records(&records, AnalyzeOptions::default()).is_err());
        let options = AnalyzeOptions {
            include_pseudo: true,
            exclude_stop_codons: false,
        };
        assert_eq!(
            analyze_records(&records, options)
                .unwrap()
                .stats
                .total_codons,
            1
        );
    }

    #[test]
    fn optionally_excludes_stop_codons_from_denominator() {
        let records = vec![record(b"ATGTAA", 1, false)];
        let options = AnalyzeOptions {
            include_pseudo: false,
            exclude_stop_codons: true,
        };
        let analysis = analyze_records(&records, options).unwrap();
        assert_eq!(analysis.stats.total_codons, 1);
        assert_eq!(analysis.stats.excluded_stop_codons, 1);
        assert_eq!(analysis.stats.count("TAA"), Some(0));
    }

    #[test]
    fn maps_codons_to_the_standard_genetic_code() {
        let stats = CodonStats::default();
        let entries: Vec<_> = stats.entries().collect();
        assert_eq!(
            entries
                .iter()
                .find(|entry| entry.codon == "ATG")
                .unwrap()
                .amino_acid,
            'M'
        );
        assert_eq!(
            entries
                .iter()
                .find(|entry| entry.codon == "AAA")
                .unwrap()
                .amino_acid,
            'K'
        );
        assert_eq!(
            entries
                .iter()
                .find(|entry| entry.codon == "TAA")
                .unwrap()
                .amino_acid,
            '*'
        );
        assert_eq!(
            entries
                .iter()
                .find(|entry| entry.codon == "TTT")
                .unwrap()
                .amino_acid,
            'F'
        );
    }

    #[test]
    fn emits_all_64_codons_in_lexicographic_order() {
        let stats = CodonStats::default();
        let entries: Vec<_> = stats.entries().collect();
        assert_eq!(entries.len(), 64);
        assert_eq!(entries[0].codon, "AAA");
        assert_eq!(entries[63].codon, "TTT");
    }
}
