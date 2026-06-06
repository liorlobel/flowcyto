#!/usr/bin/env bash
#
# Build a macOS .app bundle and a drag-to-Applications .dmg for flowcyto.
#
#   ./packaging/make-macos-app.sh
#
# Output: dist/flowcyto.app  and  dist/flowcyto-<version>.dmg
#
# Notes:
#  - Builds for the host architecture (Apple Silicon here). For a universal
#    binary, add the other target (`rustup target add x86_64-apple-darwin`) and
#    set UNIVERSAL=1.
#  - The app is ad-hoc code-signed (`codesign -s -`) so it launches on Apple
#    Silicon. It is NOT notarized: on first launch, right-click → Open (or
#    System Settings → Privacy & Security → Open Anyway) to clear Gatekeeper.
set -euo pipefail

cd "$(dirname "$0")/.."          # repo root
ROOT="$(pwd)"
APP_NAME="flowcyto"
BUNDLE_ID="online.llobel.flowcyto"
ICON_PNG="packaging/icon.png"
DIST="dist"

# cargo is not on PATH in this environment; use the rustup toolchain.
TOOLCHAIN_BIN="$HOME/.rustup/toolchains/stable-aarch64-apple-darwin/bin"
if command -v cargo >/dev/null 2>&1; then
  CARGO="cargo"
elif [ -x "$TOOLCHAIN_BIN/cargo" ]; then
  export PATH="$TOOLCHAIN_BIN:$PATH"
  CARGO="cargo"
else
  CARGO="rustup run stable cargo"
fi

VERSION="$(grep -m1 '^version' Cargo.toml | sed -E 's/.*"([^"]+)".*/\1/')"
echo "▸ flowcyto $VERSION — building release binary…"
$CARGO build --release

BIN="target/release/$APP_NAME"
[ -x "$BIN" ] || { echo "✗ binary not found at $BIN"; exit 1; }

echo "▸ assembling $APP_NAME.app…"
APP="$DIST/$APP_NAME.app"
rm -rf "$APP"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"
cp "$BIN" "$APP/Contents/MacOS/$APP_NAME"
chmod +x "$APP/Contents/MacOS/$APP_NAME"

# ── Icon: PNG → .iconset → .icns ──────────────────────────────────────────
if [ -f "$ICON_PNG" ]; then
  echo "▸ building icon…"
  ICONSET="$(mktemp -d)/$APP_NAME.iconset"
  mkdir -p "$ICONSET"
  for sz in 16 32 128 256 512; do
    sips -z $sz $sz       "$ICON_PNG" --out "$ICONSET/icon_${sz}x${sz}.png"      >/dev/null
    sips -z $((sz*2)) $((sz*2)) "$ICON_PNG" --out "$ICONSET/icon_${sz}x${sz}@2x.png" >/dev/null
  done
  iconutil -c icns "$ICONSET" -o "$APP/Contents/Resources/$APP_NAME.icns"
  ICON_LINE="<key>CFBundleIconFile</key><string>$APP_NAME</string>"
else
  echo "  (no $ICON_PNG — bundling without a custom icon)"
  ICON_LINE=""
fi

# ── Info.plist (version kept in sync with Cargo.toml) ─────────────────────
cat > "$APP/Contents/Info.plist" <<PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleName</key><string>$APP_NAME</string>
  <key>CFBundleDisplayName</key><string>flowcyto</string>
  <key>CFBundleIdentifier</key><string>$BUNDLE_ID</string>
  <key>CFBundleExecutable</key><string>$APP_NAME</string>
  $ICON_LINE
  <key>CFBundleVersion</key><string>$VERSION</string>
  <key>CFBundleShortVersionString</key><string>$VERSION</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>LSMinimumSystemVersion</key><string>11.0</string>
  <key>NSHighResolutionCapable</key><true/>
</dict>
</plist>
PLIST

printf 'APPL????' > "$APP/Contents/PkgInfo"

# ── Ad-hoc code-sign so it launches on Apple Silicon ──────────────────────
if command -v codesign >/dev/null 2>&1; then
  echo "▸ ad-hoc signing…"
  codesign --force --deep --sign - "$APP" >/dev/null 2>&1 \
    && echo "  signed (ad-hoc)" || echo "  ⚠ codesign failed — app may need right-click → Open"
fi

# ── DMG (app + Applications symlink for drag install) ─────────────────────
echo "▸ building dmg…"
DMG="$DIST/$APP_NAME-$VERSION.dmg"
rm -f "$DMG"
STAGE="$(mktemp -d)"
cp -R "$APP" "$STAGE/"
ln -s /Applications "$STAGE/Applications"
hdiutil create -volname "$APP_NAME $VERSION" -srcfolder "$STAGE" \
  -ov -format UDZO "$DMG" >/dev/null
rm -rf "$STAGE"

echo ""
echo "✓ done:"
echo "    $ROOT/$APP"
echo "    $ROOT/$DMG"
echo ""
echo "  Install: open the .dmg, drag flowcyto into Applications."
echo "  First launch: right-click flowcyto → Open (unsigned/un-notarized app)."
