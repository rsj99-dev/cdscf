use crate::error::{CdscfError, Result};
use crate::location::reverse_complement;
use crate::model::{CdsFeature, Record};
use roxmltree::Document;

const DNA_BLOCK: u8 = 0x00;
const NOTES_BLOCK: u8 = 0x06;
const COOKIE_BLOCK: u8 = 0x09;
const FEATURES_BLOCK: u8 = 0x0a;

pub(crate) fn looks_like(bytes: &[u8]) -> bool {
    if bytes.first() == Some(&COOKIE_BLOCK) {
        return true;
    }
    if bytes.len() < 6 {
        return false;
    }
    let block_type = bytes[0];
    let length = u32::from_be_bytes([bytes[1], bytes[2], bytes[3], bytes[4]]) as usize;
    let Some(end) = length.checked_add(5).filter(|end| *end <= bytes.len()) else {
        return false;
    };
    match block_type {
        DNA_BLOCK => length > 1,
        NOTES_BLOCK | FEATURES_BLOCK => {
            bytes[5..end]
                .iter()
                .copied()
                .find(|byte| !byte.is_ascii_whitespace())
                == Some(b'<')
        }
        _ => false,
    }
}

pub(crate) fn parse(bytes: &[u8]) -> Result<Vec<Record>> {
    let mut offset = 0usize;
    let mut id = None;
    let mut sequence: Option<Vec<u8>> = None;
    let mut circular = false;
    let mut raw_features = Vec::new();

    while offset < bytes.len() {
        if bytes.len() - offset < 5 {
            return Err(CdscfError::new("SnapGene 文件末尾存在不完整的块头"));
        }
        let block_type = bytes[offset];
        let length = u32::from_be_bytes([
            bytes[offset + 1],
            bytes[offset + 2],
            bytes[offset + 3],
            bytes[offset + 4],
        ]) as usize;
        offset += 5;
        let end = offset
            .checked_add(length)
            .filter(|end| *end <= bytes.len())
            .ok_or_else(|| CdscfError::new("SnapGene 块长度超出文件大小"))?;
        let data = &bytes[offset..end];
        offset = end;

        match block_type {
            NOTES_BLOCK => {
                if let Ok(xml) = std::str::from_utf8(data) {
                    id = parse_notes_name(xml);
                }
            }
            DNA_BLOCK => {
                if data.is_empty() {
                    return Err(CdscfError::new("SnapGene DNA 块为空"));
                }
                circular = data[0] & 0x01 != 0;
                sequence = Some(
                    data[1..]
                        .iter()
                        .filter(|byte| byte.is_ascii_alphabetic())
                        .map(|byte| match byte.to_ascii_uppercase() {
                            b'U' => b'T',
                            base => base,
                        })
                        .collect(),
                );
            }
            FEATURES_BLOCK => {
                let xml = std::str::from_utf8(data)
                    .map_err(|_| CdscfError::new("SnapGene 特征块不是有效 UTF-8 XML"))?;
                raw_features.extend(parse_features(xml)?);
            }
            _ => {}
        }
    }

    let sequence = sequence.ok_or_else(|| CdscfError::new("SnapGene 文件没有 DNA 序列块"))?;
    if sequence.is_empty() {
        return Err(CdscfError::new("SnapGene DNA 序列为空"));
    }
    let mut cds_features = Vec::new();
    for raw in raw_features {
        let reverse_order = raw.segments.first().is_some_and(|segment| segment.reverse);
        let mut parts = Vec::new();
        for segment in &raw.segments {
            let mut part = extract_segment(&sequence, segment, circular)?;
            if segment.reverse {
                reverse_complement(&mut part);
            }
            parts.push(part);
        }
        if reverse_order {
            parts.reverse();
        }
        let feature_sequence = parts.concat();
        cds_features.push(CdsFeature {
            name: raw.name,
            location: raw.location,
            sequence: feature_sequence,
            codon_start: raw.codon_start,
            is_pseudo: raw.is_pseudo,
        });
    }

    Ok(vec![Record {
        id: id.unwrap_or_else(|| "SnapGene".to_owned()),
        sequence,
        cds_features,
    }])
}

fn parse_notes_name(xml: &str) -> Option<String> {
    let document = Document::parse(xml).ok()?;
    document
        .descendants()
        .find(|node| node.tag_name().name() == "Name")
        .and_then(|node| node.text())
        .map(str::trim)
        .filter(|name| !name.is_empty())
        .map(str::to_owned)
}

struct RawCds {
    name: String,
    location: String,
    segments: Vec<RawSegment>,
    codon_start: usize,
    is_pseudo: bool,
}

struct RawSegment {
    start: usize,
    end: usize,
    reverse: bool,
}

fn parse_features(xml: &str) -> Result<Vec<RawCds>> {
    let document = Document::parse(xml)
        .map_err(|error| CdscfError::new(format!("SnapGene 特征 XML 无效：{error}")))?;
    let mut result = Vec::new();
    for feature in document
        .descendants()
        .filter(|node| node.tag_name().name() == "Feature")
    {
        if feature.attribute("type") != Some("CDS") {
            continue;
        }
        let feature_direction = feature.attribute("directionality").unwrap_or("1");
        let name = feature.attribute("name").unwrap_or("CDS").to_owned();
        let mut segments = Vec::new();
        let mut ranges = Vec::new();
        for segment in feature
            .children()
            .filter(|node| node.tag_name().name() == "Segment")
        {
            let range = segment.attribute("range").ok_or_else(|| {
                CdscfError::new(format!("SnapGene CDS“{name}”存在无 range 的 Segment"))
            })?;
            let (start, end) = parse_range(range, &name)?;
            let direction = segment
                .attribute("directionality")
                .unwrap_or(feature_direction);
            let reverse = match direction {
                "1" | "0" | "3" => false,
                "2" => true,
                other => {
                    return Err(CdscfError::new(format!(
                        "SnapGene CDS“{name}”具有未知 directionality={other}"
                    )))
                }
            };
            segments.push(RawSegment {
                start,
                end,
                reverse,
            });
            ranges.push(range.to_owned());
        }
        if segments.is_empty() {
            return Err(CdscfError::new(format!(
                "SnapGene CDS“{name}”不包含 Segment"
            )));
        }
        let codon_start = feature
            .attribute("codon_start")
            .or_else(|| qualifier_value(feature, "codon_start"))
            .map(|value| {
                value
                    .parse::<usize>()
                    .ok()
                    .filter(|offset| (1..=3).contains(offset))
                    .ok_or_else(|| {
                        CdscfError::new(format!("SnapGene CDS“{name}”具有无效 codon_start={value}"))
                    })
            })
            .transpose()?
            .unwrap_or(1);
        let is_pseudo = feature.attribute("pseudo").is_some()
            || feature.attribute("pseudogene").is_some()
            || has_qualifier(feature, "pseudo")
            || has_qualifier(feature, "pseudogene");
        let location = if ranges.len() == 1 {
            ranges.remove(0)
        } else {
            format!("join({})", ranges.join(","))
        };
        result.push(RawCds {
            name,
            location,
            segments,
            codon_start,
            is_pseudo,
        });
    }
    Ok(result)
}

fn qualifier_value<'a>(feature: roxmltree::Node<'a, '_>, name: &str) -> Option<&'a str> {
    feature
        .children()
        .find(|node| node.tag_name().name() == "Q" && node.attribute("name") == Some(name))
        .and_then(|node| node.children().find(|child| child.tag_name().name() == "V"))
        .and_then(|node| node.attribute("text").or_else(|| node.attribute("int")))
}

fn has_qualifier(feature: roxmltree::Node<'_, '_>, name: &str) -> bool {
    feature
        .children()
        .any(|node| node.tag_name().name() == "Q" && node.attribute("name") == Some(name))
}

fn parse_range(range: &str, name: &str) -> Result<(usize, usize)> {
    let (start, end) = range
        .split_once('-')
        .ok_or_else(|| CdscfError::new(format!("SnapGene CDS“{name}”区间“{range}”无效")))?;
    let start = start.parse::<usize>().ok().filter(|value| *value > 0);
    let end = end.parse::<usize>().ok().filter(|value| *value > 0);
    match (start, end) {
        (Some(start), Some(end)) => Ok((start, end)),
        _ => Err(CdscfError::new(format!(
            "SnapGene CDS“{name}”区间“{range}”必须包含正整数"
        ))),
    }
}

fn extract_segment(sequence: &[u8], segment: &RawSegment, circular: bool) -> Result<Vec<u8>> {
    if segment.start <= segment.end {
        if segment.end > sequence.len() {
            return Err(CdscfError::new(format!(
                "SnapGene CDS 区间 {}-{} 超出长度为 {} 的序列范围",
                segment.start,
                segment.end,
                sequence.len()
            )));
        }
        return Ok(sequence[segment.start - 1..segment.end].to_vec());
    }
    if circular && segment.start <= sequence.len() && segment.end <= sequence.len() {
        let mut result = sequence[segment.start - 1..].to_vec();
        result.extend(&sequence[..segment.end]);
        return Ok(result);
    }
    Err(CdscfError::new(format!(
        "SnapGene CDS 跨原点区间 {}-{} 无效；仅环状序列支持 start > end",
        segment.start, segment.end
    )))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn packet(block_type: u8, data: &[u8]) -> Vec<u8> {
        let mut result = vec![block_type];
        result.extend((data.len() as u32).to_be_bytes());
        result.extend(data);
        result
    }

    fn snapgene(sequence: &[u8], features: &str, circular: bool) -> Vec<u8> {
        let mut bytes = packet(COOKIE_BLOCK, b"SnapGene");
        let mut dna = vec![u8::from(circular)];
        dna.extend(sequence);
        bytes.extend(packet(DNA_BLOCK, &dna));
        bytes.extend(packet(FEATURES_BLOCK, features.as_bytes()));
        bytes
    }

    #[test]
    fn reads_forward_cds_and_ignores_other_feature_types() {
        let xml = r#"<Features><Feature type="promoter" name="noncoding"><Segment range="1-3"/></Feature><Feature type="CDS" name="gene"><Segment range="4-12"/></Feature></Features>"#;
        let bytes = snapgene(b"AAATGAAATAAGGG", xml, false);
        let records = parse(&bytes).unwrap();
        assert_eq!(records[0].cds_features.len(), 1);
        assert_eq!(records[0].cds_features[0].sequence, b"TGAAATAAG");
    }

    #[test]
    fn reads_reverse_and_joined_segments() {
        let xml = r#"<Features><Feature type="CDS" directionality="2"><Segment range="1-3"/><Segment range="7-9"/></Feature></Features>"#;
        let bytes = snapgene(b"AAACCCGGG", xml, false);
        let records = parse(&bytes).unwrap();
        assert_eq!(records[0].cds_features[0].sequence, b"CCCTTT");
    }

    #[test]
    fn supports_circular_origin_spanning_features() {
        let xml = r#"<Features><Feature type="CDS"><Segment range="8-3"/></Feature></Features>"#;
        let bytes = snapgene(b"AAACCCGGG", xml, true);
        let records = parse(&bytes).unwrap();
        assert_eq!(records[0].cds_features[0].sequence, b"GGAAA");
    }

    #[test]
    fn reads_standard_cookie_and_dna_packet_types() {
        let mut cookie = b"SnapGene".to_vec();
        cookie.extend([0, 1, 0, 1, 0, 1]);
        let mut bytes = packet(0x09, &cookie);
        let mut dna = vec![0];
        dna.extend(b"ATGAAATAA");
        bytes.extend(packet(0x00, &dna));
        let xml = r#"<Features><Feature type="CDS"><Segment range="1-9"/></Feature></Features>"#;
        bytes.extend(packet(FEATURES_BLOCK, xml.as_bytes()));

        assert!(looks_like(&bytes));
        let records = parse(&bytes).unwrap();
        assert_eq!(records[0].sequence, b"ATGAAATAA");
        assert_eq!(records[0].cds_features[0].sequence, b"ATGAAATAA");
    }

    #[test]
    fn accepts_and_detects_cookie_less_snapgene_files() {
        let xml = r#"<Features><Feature type="CDS"><Segment range="1-9"/></Feature></Features>"#;
        let mut dna = vec![0];
        dna.extend(b"ATGAAATAA");
        let mut bytes = packet(DNA_BLOCK, &dna);
        bytes.extend(packet(FEATURES_BLOCK, xml.as_bytes()));

        assert!(looks_like(&bytes));
        let records = parse(&bytes).unwrap();
        assert_eq!(records[0].cds_features[0].sequence, b"ATGAAATAA");
    }

    #[test]
    fn ignores_nonstandard_cookie_contents() {
        let xml = r#"<Features><Feature type="CDS"><Segment range="1-9"/></Feature></Features>"#;
        let mut bytes = packet(COOKIE_BLOCK, b"nonstandard-cookie");
        let mut dna = vec![0];
        dna.extend(b"ATGAAATAA");
        bytes.extend(packet(DNA_BLOCK, &dna));
        bytes.extend(packet(FEATURES_BLOCK, xml.as_bytes()));

        let records = parse(&bytes).unwrap();
        assert_eq!(records[0].cds_features[0].sequence, b"ATGAAATAA");
    }

    #[test]
    fn rejects_truncated_blocks() {
        let mut bytes = packet(COOKIE_BLOCK, b"SnapGene");
        bytes.extend([DNA_BLOCK, 0, 0, 0]);
        assert!(parse(&bytes).is_err());
    }
}
