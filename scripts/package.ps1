# package.ps1 - assemble a StimLab release (Phase E, WP-M).
#
# Produces dist/StimLab-<version>-win-x64/ containing StimLab.exe (plus any runtime
# DLLs the chosen preset needs) and a README, then zips it. The default preset is
# the fully-static `windows-static`, whose StimLab.exe is self-contained (no DLLs and
# no VC++ redist required); on first launch the app self-provisions its docking
# engine + receptors into %APPDATA%/StimLab (see the Docking panel's Provision button).
#
# Usage:
#   .\scripts\package.ps1                       # build windows-static, package, zip
#   .\scripts\package.ps1 -Preset windows       # package the dynamic build (bundles DLLs)
#   .\scripts\package.ps1 -NoBuild              # package an already-built tree
param(
    [string]$Preset  = "windows-static",
    [string]$Version = "0.1.0",
    [switch]$NoBuild
)
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root

if (-not $NoBuild) {
    Write-Host "[package] building preset '$Preset'..." -ForegroundColor Cyan
    & "$root\scripts\dev-build.ps1" $Preset
    if ($LASTEXITCODE -ne 0) { throw "build failed (exit $LASTEXITCODE)" }
}

$binDir = Join-Path $root "build\$Preset\bin"
$exe    = Join-Path $binDir "StimLab.exe"
if (-not (Test-Path $exe)) { throw "StimLab.exe not found at $exe - build the '$Preset' preset first." }

$stage = Join-Path $root "dist\StimLab-$Version-win-x64"
if (Test-Path $stage) { Remove-Item -Recurse -Force $stage }
New-Item -ItemType Directory -Force -Path $stage | Out-Null

# The exe, and any runtime DLLs sitting next to it (none for windows-static; fmt/
# spdlog/sqlite3 for the dynamic presets, copied here by vcpkg's applocal step).
Copy-Item $exe $stage
$dlls = Get-ChildItem $binDir -Filter *.dll -ErrorAction SilentlyContinue
foreach ($d in $dlls) { Copy-Item $d.FullName $stage }

$readme = @"
StimLab $Version - native CNS-stimulant computational pharmacology suite (x64 Windows)

WHAT IT IS
  Predicts what a CNS-stimulant compound IS and DOES: structure & physicochemical
  properties, molecular stability, absorption/PK, ADMET/metabolism, target binding
  affinity (real AutoDock Vina docking), similarity to known substances, legal-analog
  scoring, re-runnable prep->dock workflows, and an AI assistant that drives the app.
  Analysis only - it does NOT and will not produce synthesis routes, reaction
  conditions, precursors, or manufacturability guidance.

RUNNING
  Double-click StimLab.exe. No install, no Docker, no Python, no localhost server.
  All state lives under %APPDATA%\StimLab (database, artifacts, runtime, logs).

FIRST RUN / SELF-PROVISIONING (optional, needs internet)
  Real docking needs the AutoDock Vina engine + prepared receptors. Open the Docking
  panel and click "Provision engine + receptors": StimLab downloads vina.exe
  (size-verified) and fetches + prepares the DAT/NET/SERT/TAAR1 receptors from RCSB
  into %APPDATA%\StimLab\runtime. Until then, docking shows a clearly-labeled
  descriptor estimate. The runtime self-heals: a corrupt component is re-provisioned
  (Settings > Verify + heal runtime; manifest.json is the source of truth).

AI ASSISTANT (optional)
  The offline assistant navigates and explains. For the live provider + web tools,
  use the science build and add an Anthropic API key in Settings (stored encrypted
  via Windows DPAPI).

This build: preset '$Preset'$(if ($dlls) { " (bundled DLLs: " + (($dlls | ForEach-Object { $_.Name }) -join ', ') + ")" } else { " (fully static - single self-contained exe)" }).
"@
Set-Content -Path (Join-Path $stage "README.txt") -Value $readme -Encoding ASCII

$zip = Join-Path $root "dist\StimLab-$Version-win-x64.zip"
if (Test-Path $zip) { Remove-Item -Force $zip }
Compress-Archive -Path (Join-Path $stage '*') -DestinationPath $zip

$exeSize = "{0:N1} MB" -f ((Get-Item $exe).Length / 1MB)
Write-Host "[package] staged: $stage" -ForegroundColor Green
Write-Host "[package] exe: $exeSize   bundled DLLs: $($dlls.Count)"
Write-Host "[package] zip: $zip" -ForegroundColor Green
