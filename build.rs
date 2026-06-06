// Embed the application icon into the Windows .exe (so Explorer, the taskbar,
// and shortcuts show it). No-op on every other platform.
fn main() {
    #[cfg(windows)]
    {
        let mut res = winresource::WindowsResource::new();
        res.set_icon("packaging/icon.ico");
        let _ = res.compile();
    }
}
