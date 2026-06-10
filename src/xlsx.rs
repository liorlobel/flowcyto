//! Minimal `.xlsx` export via `rust_xlsxwriter` (pure Rust — no C or runtime deps, so it
//! keeps the single-native-binary property). Pure like `report.rs`: builds the workbook
//! bytes; the caller writes them to a file.

use rust_xlsxwriter::{Format, Workbook};

/// A cell value — kept typed so numeric columns stay numbers in Excel (sortable,
/// chartable) rather than text. Non-finite numbers are written as the text `NA`.
pub enum Cell {
    Num(f64),
    Text(String),
}

/// Serialize one worksheet (a bold header row + data rows) to in-memory `.xlsx` bytes.
pub fn sheet_bytes(sheet: &str, headers: &[&str], rows: &[Vec<Cell>]) -> Result<Vec<u8>, String> {
    let mut wb = Workbook::new();
    let ws = wb.add_worksheet();
    ws.set_name(sheet).map_err(|e| e.to_string())?;
    let bold = Format::new().set_bold();
    for (c, h) in headers.iter().enumerate() {
        ws.write_string_with_format(0, c as u16, *h, &bold).map_err(|e| e.to_string())?;
    }
    for (r, row) in rows.iter().enumerate() {
        let rr = (r + 1) as u32;
        for (c, cell) in row.iter().enumerate() {
            let cc = c as u16;
            match cell {
                Cell::Num(v) if v.is_finite() => { ws.write_number(rr, cc, *v).map_err(|e| e.to_string())?; }
                Cell::Num(_) => { ws.write_string(rr, cc, "NA").map_err(|e| e.to_string())?; }
                Cell::Text(s) => { ws.write_string(rr, cc, s).map_err(|e| e.to_string())?; }
            }
        }
    }
    ws.autofit();
    wb.save_to_buffer().map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sheet_bytes_is_a_valid_xlsx() {
        let bytes = sheet_bytes(
            "Stats",
            &["population", "count", "mfi"],
            &[
                vec![Cell::Text("CD3+".into()), Cell::Num(120.0), Cell::Num(1234.5)],
                vec![Cell::Text("empty".into()), Cell::Num(0.0), Cell::Num(f64::NAN)],
            ],
        )
        .expect("writes");
        assert!(bytes.len() > 100);
        assert_eq!(&bytes[0..2], b"PK", "xlsx is a zip (PK magic)");
    }
}
