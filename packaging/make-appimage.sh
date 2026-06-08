#!/usr/bin/env bash
#
# Build a Linux AppImage for flowcyto — a universal "download, chmod +x, run"
# binary that works across distributions.
#
#   ./packaging/make-appimage.sh
#
# Output: dist/flowcyto-<version>-<arch>.AppImage
#
# Requires a built release binary, curl, and either FUSE or (set automatically
# below) APPIMAGE_EXTRACT_AND_RUN=1 for environments without FUSE (e.g. CI).
#
# Notes:
#  - flowcyto links almost nothing unusual: its GUI libraries (libGL, libxkbcommon,
#    X11/Wayland) are loaded via dlopen and provided by the host desktop, and the
#    file dialogs use the XDG desktop portal (no GTK). So linuxdeploy bundles very
#    little beyond the binary — that is expected, not a packaging error.
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT="$(pwd)"
APP=flowcyto
VERSION="$(grep -m1 '^version' Cargo.toml | sed -E 's/.*"([^"]+)".*/\1/')"
ARCH="$(uname -m)"
DIST="dist"
mkdir -p "$DIST"

BIN="target/release/$APP"
[ -x "$BIN" ] || { echo "✗ binary not found at $BIN — run: cargo build --release"; exit 1; }

echo "▸ flowcyto $VERSION ($ARCH) — assembling AppDir…"
APPDIR="$DIST/AppDir"
rm -rf "$APPDIR"
mkdir -p "$APPDIR/usr/bin" "$APPDIR/usr/share/applications"
cp "$BIN" "$APPDIR/usr/bin/$APP"
cp "packaging/linux/$APP.desktop" "$APPDIR/usr/share/applications/$APP.desktop"
for sz in 32 48 64 128 256 512; do
  d="$APPDIR/usr/share/icons/hicolor/${sz}x${sz}/apps"
  mkdir -p "$d"
  cp "packaging/linux/icons/${sz}x${sz}/$APP.png" "$d/$APP.png"
done

# linuxdeploy bundles non-system shared libraries, writes AppRun, and (with its
# appimage plugin) packs the AppDir into a single .AppImage. Fetch both tools.
fetch() { # url dest
  [ -x "$2" ] || { echo "▸ fetching $(basename "$2")…"; curl -fsSL -o "$2" "$1"; chmod +x "$2"; }
}
BASE="https://github.com/linuxdeploy"
fetch "$BASE/linuxdeploy/releases/download/continuous/linuxdeploy-$ARCH.AppImage" "$DIST/linuxdeploy-$ARCH.AppImage"
fetch "$BASE/linuxdeploy-plugin-appimage/releases/download/continuous/linuxdeploy-plugin-appimage-$ARCH.AppImage" "$DIST/linuxdeploy-plugin-appimage-$ARCH.AppImage"

# Self-extract instead of FUSE-mount (CI runners have no FUSE); applies to the
# tools here and is baked into the produced AppImage's runtime too.
export APPIMAGE_EXTRACT_AND_RUN=1
export PATH="$ROOT/$DIST:$PATH"            # so --output appimage finds the plugin
export OUTPUT="$APP-$VERSION-$ARCH.AppImage"

echo "▸ bundling dependencies + building AppImage…"
"$DIST/linuxdeploy-$ARCH.AppImage" \
  --appdir "$APPDIR" \
  --executable "$APPDIR/usr/bin/$APP" \
  --desktop-file "$APPDIR/usr/share/applications/$APP.desktop" \
  --icon-file "packaging/linux/icons/256x256/$APP.png" \
  --output appimage

mv "$OUTPUT" "$DIST/"
echo ""
echo "✓ done: $ROOT/$DIST/$OUTPUT"
echo "  Run it with:  chmod +x $OUTPUT && ./$OUTPUT"
