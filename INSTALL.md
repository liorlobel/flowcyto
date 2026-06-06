# Installing flowcyto (macOS)

flowcyto is distributed as `flowcyto-<version>.dmg`. It is signed only ad-hoc and
is **not notarized** (no paid Apple Developer account), so macOS Gatekeeper will
warn the first time you open it. This is expected — the steps below clear it once,
permanently, per machine.

## 1. Install

1. Double-click **`flowcyto-<version>.dmg`**.
2. Drag the **flowcyto** icon onto the **Applications** folder.
3. Eject the disk image.

## 2. First launch (clear Gatekeeper — pick ONE)

**A. Terminal (most reliable, works on every macOS version)**

```bash
xattr -dr com.apple.quarantine /Applications/flowcyto.app
```

Then open flowcyto normally from Applications / Launchpad. (This removes the
"downloaded from the internet" quarantine flag; you only do it once.)

**B. Right-click → Open** (macOS 14 Sonoma and earlier)

Control-click (or right-click) **flowcyto** in Applications → **Open** →
**Open** in the dialog. After that it launches normally.

**C. "Open Anyway"** (macOS 15 Sequoia and later)

Double-click flowcyto. When macOS blocks it, go to
**System Settings → Privacy & Security**, scroll to the message about flowcyto,
and click **Open Anyway**, then confirm. After that it launches normally.

## Why the warning?

Notarization (the thing that removes the warning entirely) requires a paid Apple
Developer Program membership. The app is otherwise a normal, self-contained
native binary — no installer daemon, no network calls, nothing runs in the
background. Removing quarantine simply tells Gatekeeper you trust this app.

## Uninstall

Drag **/Applications/flowcyto.app** to the Trash. flowcyto keeps no other files
(gates/sessions are saved only where you choose via the Save dialogs).
