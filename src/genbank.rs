use crate::error::{CdscfError, Result};
use crate::location;
use crate::model::{CdsFeature, Record};

pub(crate) fn parse(bytes: &[u8]) -> Result<Vec<Record>> {
    let text = std::str::from_utf8(bytes)
        .map_err(|_| CdscfError::new("GenBank 文件不是有效的 UTF-8/ASCII 文本"))?;
    let mut records = Vec::new();
    let mut builder: Option<RecordBuilder> = None;

    for line in text.lines() {
        let line = line.trim_end_matches('\r');
        if line == "//" {
            if let Some(builder) = builder.take() {
                records.push(builder.finish()?);
            }
            continue;
        }
        if line.starts_with("LOCUS ") {
            if let Some(previous) = builder.take() {
                records.push(previous.finish()?);
            }
            let id = line.split_whitespace().nth(1).unwrap_or("unknown");
            builder = Some(RecordBuilder::new(id));
            continue;
        }
        if let Some(current) = builder.as_mut() {
            current.consume(line)?;
        }
    }
    if let Some(builder) = builder.take() {
        records.push(builder.finish()?);
    }
    if records.is_empty() {
        return Err(CdscfError::new(
            "未找到 GenBank 记录；文件应以 LOCUS 行开始",
        ));
    }
    Ok(records)
}

struct RecordBuilder {
    id: String,
    sequence: Vec<u8>,
    pending_features: Vec<RawFeature>,
    current_feature: Option<RawFeature>,
    in_features: bool,
    in_origin: bool,
}

impl RecordBuilder {
    fn new(id: &str) -> Self {
        Self {
            id: id.to_owned(),
            sequence: Vec::new(),
            pending_features: Vec::new(),
            current_feature: None,
            in_features: false,
            in_origin: false,
        }
    }

    fn consume(&mut self, line: &str) -> Result<()> {
        if line.starts_with("FEATURES") {
            self.in_features = true;
            return Ok(());
        }
        if line.starts_with("ORIGIN") {
            self.in_features = false;
            self.in_origin = true;
            return Ok(());
        }
        if self.in_origin {
            self.sequence.extend(
                line.bytes()
                    .filter(|byte| byte.is_ascii_alphabetic())
                    .map(normalize_base),
            );
            return Ok(());
        }
        if !self.in_features {
            return Ok(());
        }
        if !line.starts_with("     ") {
            self.queue_feature();
            self.in_features = false;
            return Ok(());
        }

        let is_new_feature = line
            .as_bytes()
            .get(5)
            .is_some_and(|byte| !byte.is_ascii_whitespace());
        if is_new_feature {
            self.queue_feature();
            let key = line.get(5..21).unwrap_or("").trim();
            let feature_location = line.get(21..).unwrap_or("").trim();
            self.current_feature = Some(RawFeature {
                key: key.to_owned(),
                location: feature_location.to_owned(),
                qualifiers: Vec::new(),
            });
            return Ok(());
        }

        let continuation = line.get(21..).unwrap_or("").trim();
        if continuation.is_empty() {
            return Ok(());
        }
        if let Some(feature) = self.current_feature.as_mut() {
            if continuation.starts_with('/') {
                feature.qualifiers.push(continuation.to_owned());
            } else if let Some(qualifier) = feature.qualifiers.last_mut() {
                qualifier.push(' ');
                qualifier.push_str(continuation);
            } else {
                feature.location.push_str(continuation);
            }
        }
        Ok(())
    }

    fn queue_feature(&mut self) {
        if let Some(feature) = self.current_feature.take() {
            if feature.key == "CDS" {
                self.pending_features.push(feature);
            }
        }
    }

    fn build_cds(&self, feature: RawFeature) -> Result<CdsFeature> {
        if feature.location.is_empty() {
            return Err(CdscfError::new(format!(
                "记录 {} 中存在没有位置的 CDS",
                self.id
            )));
        }
        let qualifiers = parse_qualifiers(&feature.qualifiers);
        let codon_start = qualifiers
            .iter()
            .find(|(key, _)| key == "codon_start")
            .map(|(_, value)| {
                value
                    .parse::<usize>()
                    .ok()
                    .filter(|offset| (1..=3).contains(offset))
                    .ok_or_else(|| {
                        CdscfError::new(format!(
                            "记录 {} 的 CDS“{}”具有无效 codon_start={value}",
                            self.id, feature.location
                        ))
                    })
            })
            .transpose()?
            .unwrap_or(1);
        let is_pseudo = qualifiers
            .iter()
            .any(|(key, _)| matches!(key.as_str(), "pseudo" | "pseudogene"));
        let name = ["gene", "locus_tag", "protein_id", "product"]
            .iter()
            .find_map(|wanted| {
                qualifiers
                    .iter()
                    .find(|(key, value)| key == wanted && !value.is_empty())
                    .map(|(_, value)| value.clone())
            })
            .unwrap_or_else(|| feature.location.clone());
        let parsed_location = location::parse(&feature.location).map_err(|error| {
            CdscfError::new(format!("记录 {} 的 CDS“{name}”位置无效：{error}", self.id))
        })?;
        let sequence = parsed_location.extract(&self.sequence)?;
        Ok(CdsFeature {
            name,
            location: feature.location,
            sequence,
            codon_start,
            is_pseudo,
        })
    }

    fn finish(mut self) -> Result<Record> {
        self.queue_feature();
        if self.sequence.is_empty() {
            return Err(CdscfError::new(format!(
                "记录 {} 没有 ORIGIN 序列，无法统计 CDS",
                self.id
            )));
        }
        let mut cds_features = Vec::new();
        for feature in std::mem::take(&mut self.pending_features) {
            cds_features.push(self.build_cds(feature)?);
        }
        Ok(Record {
            id: self.id,
            sequence: self.sequence,
            cds_features,
        })
    }
}

fn normalize_base(byte: u8) -> u8 {
    match byte.to_ascii_uppercase() {
        b'U' => b'T',
        base => base,
    }
}

struct RawFeature {
    key: String,
    location: String,
    qualifiers: Vec<String>,
}

fn parse_qualifiers(lines: &[String]) -> Vec<(String, String)> {
    lines
        .iter()
        .map(|line| {
            let line = line.trim_start_matches('/');
            let (key, value) = line.split_once('=').unwrap_or((line, ""));
            let value = value.trim();
            let value = if value.len() >= 2 && value.starts_with('"') && value.ends_with('"') {
                &value[1..value.len() - 1]
            } else {
                value
            };
            (key.trim().to_owned(), value.replace("\"\"", "\""))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn genbank(location: &str, qualifiers: &[&str], sequence: &str) -> Vec<u8> {
        let mut text = format!(
            "LOCUS       TEST                    {length} bp    DNA     linear\nFEATURES             Location/Qualifiers\n     CDS             {location}\n",
            length = sequence.len()
        );
        for qualifier in qualifiers {
            text.push_str(&format!("                     {qualifier}\n"));
        }
        text.push_str("ORIGIN\n        1 ");
        text.push_str(sequence);
        text.push_str("\n//\n");
        text.into_bytes()
    }

    #[test]
    fn extracts_only_the_annotated_cds() {
        let bytes = genbank("4..12", &["/gene=\"test\""], "aaaatgaaataaggg");
        let records = parse(&bytes).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "TEST");
        assert_eq!(records[0].cds_features[0].name, "test");
        assert_eq!(records[0].cds_features[0].sequence, b"ATGAAATAA");
    }

    #[test]
    fn parses_multiline_locations_and_codons_start() {
        let text = b"LOCUS       TEST                    12 bp    DNA     linear\nFEATURES             Location/Qualifiers\n     CDS             join(1..3,\n                     5..10)\n                     /codon_start=2\nORIGIN\n        1 aaacccgggttt\n//\n";
        let records = parse(text).unwrap();
        assert_eq!(records[0].cds_features[0].sequence, b"AAACCGGGT");
        assert_eq!(records[0].cds_features[0].codon_start, 2);
    }

    #[test]
    fn waits_for_origin_before_extracting_multiple_features() {
        let text = b"LOCUS       TEST                    15 bp    DNA     linear\nFEATURES             Location/Qualifiers\n     source          1..15\n     CDS             4..12\n                     /gene=\"first\"\n     gene            1..3\n                     /gene=\"noncoding\"\n     CDS             complement(1..3)\n                     /gene=\"second\"\nORIGIN\n        1 aaaatgaaataaggg\n//\n";
        let records = parse(text).unwrap();
        assert_eq!(records[0].cds_features.len(), 2);
        assert_eq!(records[0].cds_features[0].sequence, b"ATGAAATAA");
        assert_eq!(records[0].cds_features[1].sequence, b"TTT");
    }

    #[test]
    fn recognizes_pseudo_genes() {
        let bytes = genbank("1..9", &["/pseudo"], "atgaaataa");
        let records = parse(&bytes).unwrap();
        assert!(records[0].cds_features[0].is_pseudo);
    }

    #[test]
    fn parses_multiple_records() {
        let mut bytes = genbank("1..3", &[], "atg");
        bytes.extend(genbank("1..3", &[], "aaa"));
        let records = parse(&bytes).unwrap();
        assert_eq!(records.len(), 2);
    }

    #[test]
    fn rejects_remote_or_malformed_locations() {
        let bytes = genbank("OTHER.1:1..3", &[], "atg");
        assert!(parse(&bytes).is_err());
    }
}
