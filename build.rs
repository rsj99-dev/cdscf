use std::env;

fn main() {
    println!("cargo:rerun-if-changed=assets/cdscf.ico");
    println!("cargo:rerun-if-env-changed=CARGO_FEATURE_GUI");

    if env::var_os("CARGO_FEATURE_GUI").is_none()
        || env::var("CARGO_CFG_TARGET_OS").as_deref() != Ok("windows")
    {
        return;
    }

    let version = format!("{}.0", env!("CARGO_PKG_VERSION"));
    let mut resource = winresource::WindowsResource::new();
    resource
        .set_icon("assets/cdscf.ico")
        .set("FileDescription", "cdscf - CDS Codon Frequency")
        .set("ProductName", "cdscf")
        .set("InternalName", "cdscf-gui")
        .set("OriginalFilename", "cdscf-gui.exe")
        .set("FileVersion", &version)
        .set("ProductVersion", &version);

    if let Err(error) = resource.compile() {
        println!("cargo:warning=Windows 资源编译器不可用，EXE 图标未在构建期嵌入：{error}");
    }
}
