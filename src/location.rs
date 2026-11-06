use crate::error::{CdscfError, Result};

#[derive(Clone, Debug)]
pub(crate) enum Location {
    Span { start: usize, end: usize },
    Join(Vec<Location>),
    Complement(Box<Location>),
}

pub(crate) fn parse(input: &str) -> Result<Location> {
    let compact: String = input.chars().filter(|ch| !ch.is_whitespace()).collect();
    let mut parser = Parser {
        input: compact.as_bytes(),
        position: 0,
    };
    let location = parser.parse_expression()?;
    if parser.position != parser.input.len() {
        return Err(CdscfError::new(format!(
            "无法解析位置“{input}”，多余内容位于第 {} 字节",
            parser.position + 1
        )));
    }
    Ok(location)
}

impl Location {
    pub(crate) fn extract(&self, sequence: &[u8]) -> Result<Vec<u8>> {
        match self {
            Self::Span { start, end } => {
                if *start == 0 || end < start || *end > sequence.len() {
                    return Err(CdscfError::new(format!(
                        "CDS 位置 {start}..{end} 超出长度为 {} 的序列范围",
                        sequence.len()
                    )));
                }
                Ok(sequence[*start - 1..*end].to_vec())
            }
            Self::Join(parts) => {
                let mut result = Vec::new();
                for part in parts {
                    result.extend(part.extract(sequence)?);
                }
                Ok(result)
            }
            Self::Complement(part) => {
                let mut result = part.extract(sequence)?;
                reverse_complement(&mut result);
                Ok(result)
            }
        }
    }
}

pub(crate) fn reverse_complement(sequence: &mut [u8]) {
    sequence.reverse();
    for base in sequence {
        *base = match *base {
            b'A' => b'T',
            b'T' | b'U' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            b'N' => b'N',
            _ => b'N',
        };
    }
}

struct Parser<'a> {
    input: &'a [u8],
    position: usize,
}

impl Parser<'_> {
    fn parse_expression(&mut self) -> Result<Location> {
        if self.consume_keyword("complement") {
            self.expect(b'(')?;
            let inner = self.parse_expression()?;
            self.expect(b')')?;
            return Ok(Location::Complement(Box::new(inner)));
        }
        if self.consume_keyword("join") || self.consume_keyword("order") {
            self.expect(b'(')?;
            let mut parts = vec![self.parse_expression()?];
            while self.consume(b',') {
                parts.push(self.parse_expression()?);
            }
            self.expect(b')')?;
            return Ok(Location::Join(parts));
        }
        self.parse_span()
    }

    fn parse_span(&mut self) -> Result<Location> {
        self.consume_fuzzy_marker();
        let start = self.parse_number()?;
        let end = if self.consume(b'.') {
            self.expect(b'.')?;
            self.consume_fuzzy_marker();
            self.parse_number()?
        } else {
            start
        };
        if start == 0 || end < start {
            return Err(CdscfError::new(format!("无效的 CDS 区间 {start}..{end}")));
        }
        Ok(Location::Span { start, end })
    }

    fn parse_number(&mut self) -> Result<usize> {
        let start = self.position;
        while self
            .input
            .get(self.position)
            .is_some_and(u8::is_ascii_digit)
        {
            self.position += 1;
        }
        if start == self.position {
            return Err(CdscfError::new(format!(
                "位置中需要整数，位于第 {} 字节",
                self.position + 1
            )));
        }
        std::str::from_utf8(&self.input[start..self.position])
            .ok()
            .and_then(|value| value.parse().ok())
            .ok_or_else(|| CdscfError::new("CDS 坐标整数过大或无效"))
    }

    fn consume_keyword(&mut self, keyword: &str) -> bool {
        if self.input[self.position..].starts_with(keyword.as_bytes()) {
            self.position += keyword.len();
            true
        } else {
            false
        }
    }

    fn consume_fuzzy_marker(&mut self) {
        if self
            .input
            .get(self.position)
            .is_some_and(|byte| matches!(byte, b'<' | b'>'))
        {
            self.position += 1;
        }
    }

    fn consume(&mut self, expected: u8) -> bool {
        if self.input.get(self.position) == Some(&expected) {
            self.position += 1;
            true
        } else {
            false
        }
    }

    fn expect(&mut self, expected: u8) -> Result<()> {
        if self.consume(expected) {
            Ok(())
        } else {
            Err(CdscfError::new(format!(
                "位置表达式第 {} 字节应为“{}”",
                self.position + 1,
                char::from(expected)
            )))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extracts_simple_and_single_base_spans() {
        let location = parse("2..4").unwrap();
        assert_eq!(location.extract(b"ACGT").unwrap(), b"CGT");
        let location = parse("3").unwrap();
        assert_eq!(location.extract(b"ACGT").unwrap(), b"G");
    }

    #[test]
    fn complement_of_join_reverses_the_complete_join() {
        let location = parse("complement(join(1..3,7..9))").unwrap();
        assert_eq!(location.extract(b"AAACCCGGG").unwrap(), b"CCCTTT");
    }

    #[test]
    fn complemented_parts_keep_their_join_order() {
        let location = parse("join(complement(1..3),complement(7..9))").unwrap();
        assert_eq!(location.extract(b"AAACCCGGG").unwrap(), b"TTTCCC");
    }

    #[test]
    fn accepts_fuzzy_boundaries() {
        let location = parse("<2..>4").unwrap();
        assert_eq!(location.extract(b"ACGT").unwrap(), b"CGT");
    }

    #[test]
    fn rejects_out_of_bounds_coordinates() {
        let error = parse("2..5").unwrap().extract(b"ACGT").unwrap_err();
        assert!(error.to_string().contains("超出"));
    }
}
