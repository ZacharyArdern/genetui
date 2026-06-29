//! Compressed Kitty graphics protocol transmitter using raw RGBA + zlib.
//!
//! Two rendering paths are selected at runtime:
//!
//! **Unicode-placeholder path** (real Kitty terminal, `KITTY_WINDOW_ID` set):
//! Uses `U=1` + `\u{10EEEE}` diacritics so ratatui's diff engine clears old
//! content and there is no ghost/smear effect when rotating.
//!
//! **Direct-placement path** (xterm.js / VSCode and other terminals with basic
//! Kitty graphics support but no unicode placeholder support):
//! Transmits with `a=T` (no `U=1`) and `c=`/`r=` for scaling, placing the
//! image at the cursor position via save/restore.  Other cells are set to
//! skip so the image is not overwritten.  The `\u{10EEEE}` placeholder chars
//! that caused empty-box rendering in VSCode are avoided entirely.

use std::fmt::Write;
use std::io::Write as IoWrite;

use base64::Engine;
use flate2::Compression;
use flate2::write::ZlibEncoder;
use image::DynamicImage;
use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use ratatui::widgets::Widget;

/// Fixed image ID used for all Kitty image transmissions.  By reusing the
/// same ID every frame with `a=T,U=1`, Kitty atomically replaces the old
/// image data — no flicker, no delete-before-draw gap.
const IMAGE_ID: u32 = 1;

/// Returns true when running inside a real Kitty terminal that supports the
/// unicode placeholder mechanism.  VSCode/xterm.js and other basic Kitty
/// implementations do not set `KITTY_WINDOW_ID`.
fn is_kitty_native() -> bool {
    std::env::var("KITTY_WINDOW_ID").is_ok()
        || std::env::var("TERM").map(|t| t == "xterm-kitty").unwrap_or(false)
}

pub struct KittyPngImage {
    transmit: String,
    unique_id: u32,
    area: Rect,
    /// True → unicode placeholder render; false → direct-placement render.
    use_placeholders: bool,
}

impl KittyPngImage {
    /// Delete the Kitty image managed by this module.
    ///
    /// Call this when leaving FullHD mode to free Kitty-side resources and
    /// prevent stale images from lingering in the terminal's memory.
    /// Returns the escape sequence as a `String` that must be written to
    /// the terminal (e.g. via stdout).
    /// Build an iTerm2 inline-image sequence for post-draw placement in VSCode.
    ///
    /// VSCode's xterm.js has supported the iTerm2 OSC 1337 inline image protocol
    /// since 2021, predating and being more complete than its Kitty graphics support.
    /// The image is PNG-encoded (browser-native decoding), positioned with DEC
    /// save/restore cursor, and sized in terminal cells via `width=` / `height=`.
    ///
    /// Write the returned bytes directly to stdout after `terminal.draw()` so
    /// ratatui's frame cannot overwrite the image layer.
    /// Build post-draw image sequences for non-Kitty terminals.
    /// Returns `(iterm2_bel, iterm2_st, minimal)` — three variants to try in order.
    pub fn build_direct_seq(img: &DynamicImage, area: Rect) -> Option<Vec<u8>> {
        // Downscale to 300×200 to stay well under xterm.js's OSC buffer limit.
        let small = img.resize(300, 200, image::imageops::FilterType::Triangle);
        let mut png_bytes: Vec<u8> = Vec::new();
        small.write_to(&mut std::io::Cursor::new(&mut png_bytes), image::ImageFormat::Png).ok()?;
        let b64 = base64::engine::general_purpose::STANDARD.encode(&png_bytes);

        let cols = area.width;
        let rows = area.height;
        let row1 = area.top() + 1;
        let col1 = area.left() + 1;

        // Try three variants back-to-back so at least one hits if supported.
        // Variant A: full params + BEL terminator
        // Variant B: full params + ST terminator
        // Variant C: minimal (no params) — just inline=1 + data, at current cursor
        let mut seq = String::with_capacity(b64.len() * 3 + 256);

        // A — positioned, BEL
        write!(seq, "\x1b7\x1b[{row1};{col1}H").unwrap();
        write!(seq, "\x1b]1337;File=inline=1;width={cols};height={rows}:{b64}\x07").unwrap();
        write!(seq, "\x1b8").unwrap();

        // B — positioned, ST
        write!(seq, "\x1b7\x1b[{row1};{col1}H").unwrap();
        write!(seq, "\x1b]1337;File=inline=1;width={cols};height={rows}:{b64}\x1b\\").unwrap();
        write!(seq, "\x1b8").unwrap();

        // C — minimal, no positioning (wherever cursor is)
        write!(seq, "\x1b]1337;File=inline=1:{b64}\x07").unwrap();

        Some(seq.into_bytes())
    }

    pub fn cleanup_escape() -> String {
        format!("\x1b_Gq=2,a=d,d=I,i={IMAGE_ID}\x1b\\")
    }

    /// Create a new compressed Kitty image widget.
    ///
    /// The image is immediately converted to raw RGBA, zlib-compressed
    /// (level 1 / fastest), and base64-chunked into the Kitty escape
    /// sequence.  Call this outside the draw closure if you want to time
    /// encoding separately.
    ///
    /// Uses a single fixed image ID (`IMAGE_ID`).  Transmitting with
    /// `a=T,U=1` for the same ID causes Kitty to atomically replace the
    /// old image data, so there is never a visible gap between frames.
    /// No delete commands are emitted during normal rendering.
    pub fn new(img: &DynamicImage, area: Rect) -> Option<Self> {
        let (w, h) = (img.width(), img.height());
        let use_placeholders = is_kitty_native();

        let rgba = img.to_rgba8();
        let raw_bytes = rgba.as_raw();

        let mut encoder = ZlibEncoder::new(Vec::new(), Compression::fast());
        if encoder.write_all(raw_bytes).is_err() { return None; }
        let compressed = match encoder.finish() { Ok(b) => b, Err(_) => return None };

        let b64 = base64::engine::general_purpose::STANDARD.encode(&compressed);

        const CHARS_PER_CHUNK: usize = 4096;
        let chunks: Vec<&str> = b64
            .as_bytes()
            .chunks(CHARS_PER_CHUNK)
            .map(|c| std::str::from_utf8(c).expect("base64 is always valid UTF-8"))
            .collect();
        let chunk_count = chunks.len();

        let cols = area.width;
        let rows = area.height;

        let mut transmit = String::with_capacity(b64.len() + chunk_count * 80);
        for (i, chunk) in chunks.iter().enumerate() {
            write!(transmit, "\x1b_Gq=2,").unwrap();
            if i == 0 {
                if use_placeholders {
                    // Full Kitty: unicode placeholder mode — each cell carries
                    // row/col diacritics; Kitty assembles the image from them.
                    write!(transmit,
                        "i={IMAGE_ID},a=T,U=1,f=32,o=z,t=d,s={w},v={h},c={cols},r={rows},"
                    ).unwrap();
                } else {
                    // Basic Kitty (VSCode/xterm.js): direct placement at cursor.
                    // c=/r= tell the terminal to scale the image to that many cells.
                    // No U=1 → no unicode placeholder characters in the buffer.
                    write!(transmit,
                        "i={IMAGE_ID},a=T,f=32,o=z,t=d,s={w},v={h},c={cols},r={rows},"
                    ).unwrap();
                }
            }
            let more = u8::from(chunk_count > (i + 1));
            write!(transmit, "m={more};{chunk}\x1b\\").unwrap();
        }

        Some(Self { transmit, unique_id: IMAGE_ID, area, use_placeholders })
    }
}

impl Widget for KittyPngImage {
    fn render(self, area: Rect, buf: &mut Buffer) {
        if self.use_placeholders {
            self.render_placeholders(area, buf);
        } else {
            self.render_direct(area, buf);
        }
    }
}

impl KittyPngImage {
    /// Unicode-placeholder render (real Kitty terminal).
    /// Every cell carries row/column diacritics so ratatui's diff engine
    /// properly clears stale content between frames.
    fn render_placeholders(mut self, area: Rect, buf: &mut Buffer) {
        let [id_extra, id_r, id_g, id_b] = self.unique_id.to_be_bytes();
        let id_color = format!("\x1b[38;2;{id_r};{id_g};{id_b}m");
        let full_height = area.height.min(self.area.height).min(297);
        let full_width  = area.width.min(self.area.width);

        for y in 0..full_height {
            let mut symbol = if y == 0 {
                std::mem::take(&mut self.transmit)
            } else {
                String::new()
            };
            let width_usize = usize::from(full_width);
            symbol.reserve(id_color.len() + (width_usize * 16) + 32);

            write!(
                symbol,
                "\x1b[s{id_color}\u{10EEEE}{}{}{}",
                diacritic(y), diacritic(0), diacritic(u16::from(id_extra)),
            ).unwrap();
            symbol.extend(std::iter::repeat_n('\u{10EEEE}', width_usize.saturating_sub(1)));

            let right = area.width.saturating_sub(1);
            let down  = area.height.saturating_sub(1);
            write!(symbol, "\x1b[u\x1b[{right}C\x1b[{down}B").unwrap();

            if let Some(cell) = buf.cell_mut((area.left(), area.top() + y)) {
                cell.set_symbol(&symbol);
            }
            for x in 1..full_width {
                if let Some(cell) = buf.cell_mut((area.left() + x, area.top() + y)) {
                    cell.set_skip(true);
                }
            }
        }
    }

    /// Direct-placement render (VSCode/xterm.js and other basic Kitty terminals).
    ///
    /// The Kitty escape sequence (without `U=1`) is placed into the top-left
    /// cell wrapped in cursor save/restore.  All other cells are filled with
    /// dark-background spaces rather than skipped — `set_skip` leaves cells
    /// unwritten so the terminal shows raw white on first render.  In xterm.js,
    /// Kitty images are composited on a separate canvas layer above the text
    /// layer, so these spaces appear behind the image.
    fn render_direct(mut self, area: Rect, buf: &mut Buffer) {
        use ratatui::style::{Color, Style};

        let right = area.width.saturating_sub(1);
        let down  = area.height.saturating_sub(1);

        let mut symbol = String::new();
        write!(symbol, "\x1b[s").unwrap();
        symbol.push_str(&std::mem::take(&mut self.transmit));
        write!(symbol, "\x1b[u\x1b[{right}C\x1b[{down}B").unwrap();

        if let Some(cell) = buf.cell_mut((area.left(), area.top())) {
            cell.set_symbol(&symbol);
        }

        // Fill remaining cells with dark bg so ratatui writes them to the
        // terminal, preventing the white-background bleed-through.
        let dark = Style::default().bg(Color::Rgb(18, 18, 30));
        for y in 0..area.height {
            let start_x = if y == 0 { 1 } else { 0 };
            for x in start_x..area.width {
                if let Some(cell) = buf.cell_mut((area.left() + x, area.top() + y)) {
                    cell.set_symbol(" ");
                    cell.set_style(dark);
                }
            }
        }
    }
}

/// Diacritics table for Kitty unicode placeholders.
/// From <https://sw.kovidgoyal.net/kitty/_downloads/1792bad15b12979994cd6ecc54c967a6/rowcolumn-diacritics.txt>
/// See <https://sw.kovidgoyal.net/kitty/graphics-protocol/#unicode-placeholders>
static DIACRITICS: [char; 297] = [
    '\u{305}',
    '\u{30D}',
    '\u{30E}',
    '\u{310}',
    '\u{312}',
    '\u{33D}',
    '\u{33E}',
    '\u{33F}',
    '\u{346}',
    '\u{34A}',
    '\u{34B}',
    '\u{34C}',
    '\u{350}',
    '\u{351}',
    '\u{352}',
    '\u{357}',
    '\u{35B}',
    '\u{363}',
    '\u{364}',
    '\u{365}',
    '\u{366}',
    '\u{367}',
    '\u{368}',
    '\u{369}',
    '\u{36A}',
    '\u{36B}',
    '\u{36C}',
    '\u{36D}',
    '\u{36E}',
    '\u{36F}',
    '\u{483}',
    '\u{484}',
    '\u{485}',
    '\u{486}',
    '\u{487}',
    '\u{592}',
    '\u{593}',
    '\u{594}',
    '\u{595}',
    '\u{597}',
    '\u{598}',
    '\u{599}',
    '\u{59C}',
    '\u{59D}',
    '\u{59E}',
    '\u{59F}',
    '\u{5A0}',
    '\u{5A1}',
    '\u{5A8}',
    '\u{5A9}',
    '\u{5AB}',
    '\u{5AC}',
    '\u{5AF}',
    '\u{5C4}',
    '\u{610}',
    '\u{611}',
    '\u{612}',
    '\u{613}',
    '\u{614}',
    '\u{615}',
    '\u{616}',
    '\u{617}',
    '\u{657}',
    '\u{658}',
    '\u{659}',
    '\u{65A}',
    '\u{65B}',
    '\u{65D}',
    '\u{65E}',
    '\u{6D6}',
    '\u{6D7}',
    '\u{6D8}',
    '\u{6D9}',
    '\u{6DA}',
    '\u{6DB}',
    '\u{6DC}',
    '\u{6DF}',
    '\u{6E0}',
    '\u{6E1}',
    '\u{6E2}',
    '\u{6E4}',
    '\u{6E7}',
    '\u{6E8}',
    '\u{6EB}',
    '\u{6EC}',
    '\u{730}',
    '\u{732}',
    '\u{733}',
    '\u{735}',
    '\u{736}',
    '\u{73A}',
    '\u{73D}',
    '\u{73F}',
    '\u{740}',
    '\u{741}',
    '\u{743}',
    '\u{745}',
    '\u{747}',
    '\u{749}',
    '\u{74A}',
    '\u{7EB}',
    '\u{7EC}',
    '\u{7ED}',
    '\u{7EE}',
    '\u{7EF}',
    '\u{7F0}',
    '\u{7F1}',
    '\u{7F3}',
    '\u{816}',
    '\u{817}',
    '\u{818}',
    '\u{819}',
    '\u{81B}',
    '\u{81C}',
    '\u{81D}',
    '\u{81E}',
    '\u{81F}',
    '\u{820}',
    '\u{821}',
    '\u{822}',
    '\u{823}',
    '\u{825}',
    '\u{826}',
    '\u{827}',
    '\u{829}',
    '\u{82A}',
    '\u{82B}',
    '\u{82C}',
    '\u{82D}',
    '\u{951}',
    '\u{953}',
    '\u{954}',
    '\u{F82}',
    '\u{F83}',
    '\u{F86}',
    '\u{F87}',
    '\u{135D}',
    '\u{135E}',
    '\u{135F}',
    '\u{17DD}',
    '\u{193A}',
    '\u{1A17}',
    '\u{1A75}',
    '\u{1A76}',
    '\u{1A77}',
    '\u{1A78}',
    '\u{1A79}',
    '\u{1A7A}',
    '\u{1A7B}',
    '\u{1A7C}',
    '\u{1B6B}',
    '\u{1B6D}',
    '\u{1B6E}',
    '\u{1B6F}',
    '\u{1B70}',
    '\u{1B71}',
    '\u{1B72}',
    '\u{1B73}',
    '\u{1CD0}',
    '\u{1CD1}',
    '\u{1CD2}',
    '\u{1CDA}',
    '\u{1CDB}',
    '\u{1CE0}',
    '\u{1DC0}',
    '\u{1DC1}',
    '\u{1DC3}',
    '\u{1DC4}',
    '\u{1DC5}',
    '\u{1DC6}',
    '\u{1DC7}',
    '\u{1DC8}',
    '\u{1DC9}',
    '\u{1DCB}',
    '\u{1DCC}',
    '\u{1DD1}',
    '\u{1DD2}',
    '\u{1DD3}',
    '\u{1DD4}',
    '\u{1DD5}',
    '\u{1DD6}',
    '\u{1DD7}',
    '\u{1DD8}',
    '\u{1DD9}',
    '\u{1DDA}',
    '\u{1DDB}',
    '\u{1DDC}',
    '\u{1DDD}',
    '\u{1DDE}',
    '\u{1DDF}',
    '\u{1DE0}',
    '\u{1DE1}',
    '\u{1DE2}',
    '\u{1DE3}',
    '\u{1DE4}',
    '\u{1DE5}',
    '\u{1DE6}',
    '\u{1DFE}',
    '\u{20D0}',
    '\u{20D1}',
    '\u{20D4}',
    '\u{20D5}',
    '\u{20D6}',
    '\u{20D7}',
    '\u{20DB}',
    '\u{20DC}',
    '\u{20E1}',
    '\u{20E7}',
    '\u{20E9}',
    '\u{20F0}',
    '\u{2CEF}',
    '\u{2CF0}',
    '\u{2CF1}',
    '\u{2DE0}',
    '\u{2DE1}',
    '\u{2DE2}',
    '\u{2DE3}',
    '\u{2DE4}',
    '\u{2DE5}',
    '\u{2DE6}',
    '\u{2DE7}',
    '\u{2DE8}',
    '\u{2DE9}',
    '\u{2DEA}',
    '\u{2DEB}',
    '\u{2DEC}',
    '\u{2DED}',
    '\u{2DEE}',
    '\u{2DEF}',
    '\u{2DF0}',
    '\u{2DF1}',
    '\u{2DF2}',
    '\u{2DF3}',
    '\u{2DF4}',
    '\u{2DF5}',
    '\u{2DF6}',
    '\u{2DF7}',
    '\u{2DF8}',
    '\u{2DF9}',
    '\u{2DFA}',
    '\u{2DFB}',
    '\u{2DFC}',
    '\u{2DFD}',
    '\u{2DFE}',
    '\u{2DFF}',
    '\u{A66F}',
    '\u{A67C}',
    '\u{A67D}',
    '\u{A6F0}',
    '\u{A6F1}',
    '\u{A8E0}',
    '\u{A8E1}',
    '\u{A8E2}',
    '\u{A8E3}',
    '\u{A8E4}',
    '\u{A8E5}',
    '\u{A8E6}',
    '\u{A8E7}',
    '\u{A8E8}',
    '\u{A8E9}',
    '\u{A8EA}',
    '\u{A8EB}',
    '\u{A8EC}',
    '\u{A8ED}',
    '\u{A8EE}',
    '\u{A8EF}',
    '\u{A8F0}',
    '\u{A8F1}',
    '\u{AAB0}',
    '\u{AAB2}',
    '\u{AAB3}',
    '\u{AAB7}',
    '\u{AAB8}',
    '\u{AABE}',
    '\u{AABF}',
    '\u{AAC1}',
    '\u{FE20}',
    '\u{FE21}',
    '\u{FE22}',
    '\u{FE23}',
    '\u{FE24}',
    '\u{FE25}',
    '\u{FE26}',
    '\u{10A0F}',
    '\u{10A38}',
    '\u{1D185}',
    '\u{1D186}',
    '\u{1D187}',
    '\u{1D188}',
    '\u{1D189}',
    '\u{1D1AA}',
    '\u{1D1AB}',
    '\u{1D1AC}',
    '\u{1D1AD}',
    '\u{1D242}',
    '\u{1D243}',
    '\u{1D244}',
];

#[inline]
fn diacritic(y: u16) -> char {
    *DIACRITICS.get(usize::from(y)).unwrap_or(&DIACRITICS[0])
}
