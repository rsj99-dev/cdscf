use crate::error::Result;
use crate::stats::Analysis;
use std::io::Write;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum OutputFormat {
    #[default]
    Table,
    Csv,
    Json,
}

pub fn write_analysis(
    mut writer: impl Write,
    analysis: &Analysis,
    format: OutputFormat,
) -> Result<()> {
    match format {
        OutputFormat::Table => write_table(&mut writer, analysis)?,
        OutputFormat::Csv => write_csv(&mut writer, analysis)?,
        OutputFormat::Json => write_json(&mut writer, analysis)?,
    }
    writer.flush()?;
    Ok(())
}

fn write_table(mut writer: impl Write, analysis: &Analysis) -> Result<()> {
    writeln!(writer, "记录数: {}", analysis.record_count())?;
    writeln!(writer, "记录 ID: {}", analysis.record_ids.join(", "))?;
    writeln!(writer, "统计的 CDS 数: {}", analysis.cds_features)?;
    writeln!(writer, "跳过的伪基因 CDS 数: {}", analysis.skipped_pseudo)?;
    writeln!(writer, "有效密码子总数: {}", analysis.stats.total_codons)?;
    writeln!(
        writer,
        "跳过的模糊密码子数: {}",
        analysis.stats.ambiguous_codons
    )?;
    writeln!(
        writer,
        "不足三联体的末端碱基数: {}",
        analysis.stats.trailing_bases
    )?;
    writeln!(
        writer,
        "排除的终止密码子数: {}",
        analysis.stats.excluded_stop_codons
    )?;
    writeln!(writer)?;
    writeln!(
        writer,
        "{:<6} {:<4} {:>12} {:>16}",
        "CODON", "AA", "COUNT", "FREQUENCY"
    )?;
    for entry in analysis.stats.entries() {
        writeln!(
            writer,
            "{:<6} {:<4} {:>12} {:>16.8}",
            entry.codon, entry.amino_acid, entry.count, entry.frequency
        )?;
    }
    Ok(())
}

fn write_csv(mut writer: impl Write, analysis: &Analysis) -> Result<()> {
    writeln!(writer, "codon,amino_acid,count,frequency")?;
    for entry in analysis.stats.entries() {
        writeln!(
            writer,
            "{},{},{},{:.8}",
            entry.codon, entry.amino_acid, entry.count, entry.frequency
        )?;
    }
    Ok(())
}

fn write_json(mut writer: impl Write, analysis: &Analysis) -> Result<()> {
    writeln!(writer, "{{")?;
    write!(writer, "  \"records\": [")?;
    for (index, id) in analysis.record_ids.iter().enumerate() {
        if index > 0 {
            write!(writer, ", ")?;
        }
        write!(writer, "\"{}\"", json_escape(id))?;
    }
    writeln!(writer, "],")?;
    writeln!(writer, "  \"cds_features\": {},", analysis.cds_features)?;
    writeln!(writer, "  \"skipped_pseudo\": {},", analysis.skipped_pseudo)?;
    writeln!(
        writer,
        "  \"total_codons\": {},",
        analysis.stats.total_codons
    )?;
    writeln!(
        writer,
        "  \"ambiguous_codons\": {},",
        analysis.stats.ambiguous_codons
    )?;
    writeln!(
        writer,
        "  \"trailing_bases\": {},",
        analysis.stats.trailing_bases
    )?;
    writeln!(
        writer,
        "  \"excluded_stop_codons\": {},",
        analysis.stats.excluded_stop_codons
    )?;
    writeln!(writer, "  \"codons\": [")?;
    for (index, entry) in analysis.stats.entries().enumerate() {
        writeln!(
            writer,
            "    {{\"codon\": \"{}\", \"amino_acid\": \"{}\", \"count\": {}, \"frequency\": {:.8}}}{}",
            entry.codon,
            entry.amino_acid,
            entry.count,
            entry.frequency,
            if index == 63 { "" } else { "," }
        )?;
    }
    writeln!(writer, "  ]")?;
    writeln!(writer, "}}")?;
    Ok(())
}

fn json_escape(value: &str) -> String {
    let mut escaped = String::new();
    for character in value.chars() {
        match character {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            ch if ch.is_control() => {
                use std::fmt::Write as _;
                let _ = write!(escaped, "\\u{:04x}", ch as u32);
            }
            ch => escaped.push(ch),
        }
    }
    escaped
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::CodonStats;

    #[test]
    fn csv_has_header_and_64_rows() {
        let analysis = Analysis {
            record_ids: vec!["test".to_owned()],
            cds_features: 0,
            skipped_pseudo: 0,
            stats: CodonStats::default(),
        };
        let mut output = Vec::new();
        write_analysis(&mut output, &analysis, OutputFormat::Csv).unwrap();
        assert_eq!(String::from_utf8(output).unwrap().lines().count(), 65);
    }

    #[test]
    fn json_escapes_record_ids() {
        assert_eq!(json_escape("a\"b\nc"), "a\\\"b\\nc");
    }
}
