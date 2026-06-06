; Inno Setup script for flowcyto (Windows installer).
; Built in CI on a Windows runner:
;   ISCC.exe /DMyAppVersion=<x.y.z> packaging\windows\flowcyto.iss
; Produces packaging\windows\Output\flowcyto-<version>-setup.exe

#define MyAppName "flowcyto"
#ifndef MyAppVersion
  #define MyAppVersion "0.0.0"
#endif
#define MyAppExe "flowcyto.exe"
#define MyAppPublisher "Lior Lobel"
#define MyAppURL "https://github.com/liorlobel/flowcyto"

[Setup]
AppId={{B9D6F1C2-7A3E-4B1A-9E2C-FLOWCYTO0001}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
AppPublisherURL={#MyAppURL}
AppSupportURL={#MyAppURL}
DefaultDirName={autopf}\flowcyto
DefaultGroupName=flowcyto
DisableProgramGroupPage=yes
UninstallDisplayIcon={app}\{#MyAppExe}
SetupIconFile=..\icon.ico
OutputDir=Output
OutputBaseFilename=flowcyto-{#MyAppVersion}-setup
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Additional icons:"; Flags: unchecked

[Files]
Source: "..\..\target\release\{#MyAppExe}"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{group}\flowcyto"; Filename: "{app}\{#MyAppExe}"
Name: "{group}\Uninstall flowcyto"; Filename: "{uninstallexe}"
Name: "{autodesktop}\flowcyto"; Filename: "{app}\{#MyAppExe}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#MyAppExe}"; Description: "Launch flowcyto"; Flags: nowait postinstall skipifsilent
