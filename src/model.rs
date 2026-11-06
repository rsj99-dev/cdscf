/// 一个 GenBank 记录或一个 SnapGene 文档。
#[derive(Clone, Debug)]
pub struct Record {
    pub id: String,
    pub sequence: Vec<u8>,
    pub cds_features: Vec<CdsFeature>,
}

/// 已按 5'→3' 编码方向提取的 CDS。
#[derive(Clone, Debug)]
pub struct CdsFeature {
    pub name: String,
    pub location: String,
    pub sequence: Vec<u8>,
    pub codon_start: usize,
    pub is_pseudo: bool,
}
