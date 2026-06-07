//! Lightweight "check for updates" against the GitHub Releases API.
//!
//! Privacy: this is the app's ONLY network call, and it fires *only* when the user
//! explicitly clicks "Check for Updates" — never on launch, never in the background.
//! It downloads nothing executable; it just compares versions and (on the user's
//! click) opens the release page in the browser so they install the signed installer
//! the same way they already do.

const RELEASES_API: &str = "https://api.github.com/repos/liorlobel/flowcyto/releases/latest";
/// Human-facing releases page (fallback if the API response lacks an html_url).
pub const RELEASES_PAGE: &str = "https://github.com/liorlobel/flowcyto/releases/latest";

#[derive(Clone, Debug)]
pub struct UpdateInfo {
    pub current: String,
    pub latest: String,
    /// True if `latest` is strictly newer than the running version.
    pub newer: bool,
    /// Release page to open in the browser.
    pub url: String,
}

pub fn current_version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

/// Query the GitHub Releases API for the latest published version. **Blocking** —
/// run it on a background thread so the UI never freezes.
pub fn check_latest() -> Result<UpdateInfo, String> {
    let ua = format!("flowcyto/{}", current_version());
    // 15 s global timeout so an offline / hung click fails cleanly instead of leaving
    // the background thread (and the toolbar spinner) waiting on the OS TCP timeout.
    let agent: ureq::Agent = ureq::Agent::config_builder()
        .timeout_global(Some(std::time::Duration::from_secs(15)))
        .build()
        .into();
    let mut resp = agent
        .get(RELEASES_API)
        .header("User-Agent", &ua)
        .header("Accept", "application/vnd.github+json")
        .call()
        .map_err(|e| format!("network error: {e}"))?;
    let body = resp
        .body_mut()
        .read_to_string()
        .map_err(|e| format!("read error: {e}"))?;
    let json: serde_json::Value =
        serde_json::from_str(&body).map_err(|e| format!("unexpected response: {e}"))?;
    let tag = json
        .get("tag_name")
        .and_then(|v| v.as_str())
        .ok_or("response had no tag_name")?;
    let url = json
        .get("html_url")
        .and_then(|v| v.as_str())
        .unwrap_or(RELEASES_PAGE)
        .to_string();
    let latest = tag.trim_start_matches('v').to_string();
    let current = current_version().to_string();
    let newer = is_newer(&latest, &current);
    Ok(UpdateInfo { current, latest, newer, url })
}

/// True if `latest` is a strictly higher version than `current` (semantic X.Y.Z).
/// Unparseable input is treated as "not newer" (fail safe — never nags spuriously).
pub fn is_newer(latest: &str, current: &str) -> bool {
    match (parse_ver(latest), parse_ver(current)) {
        (Some(l), Some(c)) => l > c,
        _ => false,
    }
}

/// Parse "v1.2.3" / "1.2" / "1.2.3-rc1" into a comparable (major, minor, patch).
fn parse_ver(s: &str) -> Option<(u32, u32, u32)> {
    let s = s.trim().trim_start_matches('v');
    let core = s.split(['-', '+']).next().unwrap_or(s); // drop pre-release / build suffix
    let mut it = core.split('.');
    let major = it.next()?.parse().ok()?;
    let minor = it.next().unwrap_or("0").parse().ok()?;
    let patch = it.next().unwrap_or("0").parse().ok()?;
    Some((major, minor, patch))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn newer_detection() {
        assert!(is_newer("0.1.8", "0.1.7"));
        assert!(is_newer("0.2.0", "0.1.9"));
        assert!(is_newer("1.0.0", "0.9.9"));
        assert!(!is_newer("0.1.7", "0.1.7"));
        assert!(!is_newer("0.1.6", "0.1.7"));
    }

    #[test]
    fn parses_v_prefix_and_suffix() {
        assert_eq!(parse_ver("v0.1.7"), Some((0, 1, 7)));
        assert_eq!(parse_ver("0.1.7-rc1"), Some((0, 1, 7)));
        assert_eq!(parse_ver("1.2"), Some((1, 2, 0)));
        assert_eq!(parse_ver("garbage"), None);
    }

    #[test]
    fn malformed_is_not_newer() {
        // Never nag on unparseable versions.
        assert!(!is_newer("garbage", "0.1.7"));
        assert!(!is_newer("0.1.8", "garbage"));
    }
}
