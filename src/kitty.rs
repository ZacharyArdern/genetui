/// Write the Kitty cleanup escape sequence to stdout to delete the protein image.
pub fn kitty_delete() {
    use std::io::Write;
    let seq = crate::pv::kitty_png::KittyPngImage::cleanup_escape();
    if let Ok(mut tty) = std::fs::OpenOptions::new().write(true).open("/dev/tty") {
        let _ = tty.write_all(seq.as_bytes());
    }
}
