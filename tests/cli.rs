//! End-to-end CLI tests: run the built `flowcyto` binary against a committed
//! fixture and assert on stdout. Covers arg parsing + output formatting, the
//! one layer the in-crate unit tests can't reach.

use std::process::Command;

/// Path to the binary cargo built for this test run.
const BIN: &str = env!("CARGO_BIN_EXE_flowcyto");
const FIXTURE: &str = "tests/fixtures/tiny.fcs";

fn run(args: &[&str]) -> (String, String, bool) {
    let out = Command::new(BIN)
        .args(args)
        .output()
        .expect("failed to spawn flowcyto");
    (
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
        out.status.success(),
    )
}

#[test]
fn info_reports_dimensions_and_channels() {
    let (stdout, stderr, ok) = run(&["info", FIXTURE]);
    assert!(ok, "info failed: {stderr}");
    assert!(stdout.contains("Version : FCS3.0"), "{stdout}");
    assert!(stdout.contains("Events  : 4"), "{stdout}");
    assert!(stdout.contains("Params  : 3"), "{stdout}");
    assert!(stdout.contains("FITC-A"), "{stdout}");
    assert!(stdout.contains("CD11c"), "channel label missing: {stdout}");
    assert!(stdout.contains("$SPILLOVER"), "{stdout}");
}

#[test]
fn stats_lists_every_channel() {
    let (stdout, stderr, ok) = run(&["stats", FIXTURE]);
    assert!(ok, "stats failed: {stderr}");
    for ch in ["FSC-A", "FITC-A", "PE-A"] {
        assert!(stdout.contains(ch), "stats missing {ch}: {stdout}");
    }
}

#[test]
fn spillover_shows_off_diagonal() {
    let (stdout, stderr, ok) = run(&["spillover", FIXTURE]);
    assert!(ok, "spillover failed: {stderr}");
    assert!(stdout.contains("0.2000"), "expected 20% spillover: {stdout}");
    assert!(
        stdout.contains("real compensation matrix"),
        "{stdout}"
    );
}

#[test]
fn missing_file_errors_cleanly() {
    let (_stdout, _stderr, ok) = run(&["info", "tests/fixtures/does_not_exist.fcs"]);
    assert!(!ok, "opening a missing file should fail");
}
