use std::fmt;

/// cdscf 的统一错误类型。
#[derive(Debug)]
pub struct CdscfError {
    message: String,
}

impl CdscfError {
    pub fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
        }
    }
}

impl fmt::Display for CdscfError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.message)
    }
}

impl std::error::Error for CdscfError {}

impl From<std::io::Error> for CdscfError {
    fn from(error: std::io::Error) -> Self {
        Self::new(format!("输入/输出错误：{error}"))
    }
}

pub type Result<T> = std::result::Result<T, CdscfError>;
