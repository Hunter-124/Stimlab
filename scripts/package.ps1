# package.ps1 - assemble a BioCAD release (Phase E, WP-M).
#
# Produces dist/BioCAD-<version>-win-x64/ containing BioCAD.exe (plus any runtime
# DLLs the chosen preset needs) and a README, then zips it. The default preset is
# the fully-static `windows-static`, whose BioCAD.exe is self-contained (no DLLs and
# no VC++ redist required); on first launch the app self-provisions its docking
# engine + receptors into %APPDATA%/BioCAD (see the Docking panel's Provision button).
#
# Usage:
#   .\scripts\package.ps1                       # build windows-static, package, zip
#   .\scripts\package.ps1 -Preset windows       # package the dynamic build (bundles DLLs)
#   .\scripts\package.ps1 -NoBuild              # package an already-built tree
param(
    [string]$Preset  = "windows-static",
    [string]$Version = "",
    [switch]$NoBuild
)
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root
if (-not $Version) { $Version = & "$PSScriptRoot\version.ps1" }

if (-not $NoBuild) {
    Write-Host "[package] building preset '$Preset'..." -ForegroundColor Cyan
    & "$root\scripts\dev-build.ps1" $Preset
    if ($LASTEXITCODE -ne 0) { throw "build failed (exit $LASTEXITCODE)" }
}

$binDir = Join-Path $root "build\$Preset\bin"
$exe    = Join-Path $binDir "BioCAD.exe"
if (-not (Test-Path $exe)) { throw "BioCAD.exe not found at $exe - build the '$Preset' preset first." }

$stage = Join-Path $root "dist\BioCAD-$Version-win-x64"
if (Test-Path $stage) { Remove-Item -Recurse -Force $stage }
New-Item -ItemType Directory -Force -Path $stage | Out-Null

# The exe, and any runtime DLLs sitting next to it (none for windows-static; fmt/
# spdlog/sqlite3 for the dynamic presets, copied here by vcpkg's applocal step).
Copy-Item $exe $stage
$dlls = Get-ChildItem $binDir -Filter *.dll -ErrorAction SilentlyContinue
foreach ($d in $dlls) { Copy-Item $d.FullName $stage }

# The compound/target catalog is data, not code: without assets/packs the app has
# no library and no receptor presets, so this is a required part of the payload.
$packSrc = Join-Path $root "assets\packs"
if (-not (Test-Path $packSrc)) { throw "assets\packs not found at $packSrc - the release would ship with no catalog." }
$packDst = Join-Path $stage "assets\packs"
New-Item -ItemType Directory -Force -Path $packDst | Out-Null
Copy-Item (Join-Path $packSrc "*.json") $packDst

# The science presets statically bundle libcurl, so the live Anthropic provider + web
# tools ship inside this exe; the plain presets are curl-free (offline assistant only).
$science = $Preset -like "*science*"
$aiPara = if ($science) {
@"
AI ASSISTANT (live provider + web tools INCLUDED in this build)
  This build bundles the live Anthropic provider and the keyless web_search /
  web_fetch tools (static libcurl, Windows-native Schannel TLS). The offline
  assistant works out of the box; add an Anthropic API key in Settings (stored
  encrypted via Windows DPAPI) to enable the live agent.
"@
} else {
@"
AI ASSISTANT (optional)
  The offline assistant navigates and explains. For the live provider + web tools,
  build the science release (build.ps1 -Release -Science) and add an Anthropic API
  key in Settings (stored encrypted via Windows DPAPI).
"@
}

$readme = @"
BioCAD $Version - native workstation for molecular, protein, and pharmacological analysis (x64 Windows)

WHAT IT IS
  Analyses what a compound IS and DOES: structure & physicochemical properties,
  molecular stability, absorption/PK, ADMET/metabolism, target binding affinity
  (real AutoDock Vina docking), similarity to known substances, legal-analog
  scoring, re-runnable prep->dock workflows, and an AI assistant that drives the app.
  Every derived number carries a provenance tier, so a measurement is never shown
  as if it were a prediction. Analysis only - it does NOT and will not produce
  synthesis routes, reaction conditions, precursors, or manufacturability guidance.

RUNNING
  Double-click BioCAD.exe. No install, no Docker, no Python, no localhost server.
  All state lives under %APPDATA%\BioCAD (database, artifacts, runtime, logs).

FIRST RUN / SELF-PROVISIONING (optional, needs internet)
  Real docking needs the AutoDock Vina engine + prepared receptors. Open the Docking
  panel and click "Provision engine + receptors": BioCAD downloads vina.exe
  (size-verified) and fetches + prepares the DAT/NET/SERT/TAAR1 receptors from RCSB
  into %APPDATA%\BioCAD\runtime. Until then, docking shows a clearly-labeled
  descriptor estimate. The runtime self-heals: a corrupt component is re-provisioned
  (Settings > Verify + heal runtime; manifest.json is the source of truth).

INTEGRITY (this build is intentionally UNSIGNED)
  BioCAD is an educational / personal-use project and ships without an Authenticode
  signature, so Windows SmartScreen may show "Windows protected your PC" on first run:
  click More info > Run anyway. Nothing here is signed or certified. (A maintainer with
  a code-signing certificate can sign the exe via scripts\sign.ps1; absent a cert it is
  a clean no-op.) The downloaded docking engine is size- and content-verified, and an
  optional SHA-256 pin (BIOCAD_VINA_SHA256) cryptographically verifies vina.exe.

$aiPara
This build: preset '$Preset'$(if ($dlls) { " (bundled DLLs: " + (($dlls | ForEach-Object { $_.Name }) -join ', ') + ")" } else { " (fully static - single self-contained exe)" })$(if ($science) { " - live agent + web tools bundled" }).
"@
Set-Content -Path (Join-Path $stage "README.txt") -Value $readme -Encoding ASCII

$zip = Join-Path $root "dist\BioCAD-$Version-win-x64.zip"
if (Test-Path $zip) { Remove-Item -Force $zip }
Compress-Archive -Path (Join-Path $stage '*') -DestinationPath $zip

$exeSize = "{0:N1} MB" -f ((Get-Item $exe).Length / 1MB)
Write-Host "[package] staged: $stage" -ForegroundColor Green
Write-Host "[package] exe: $exeSize   bundled DLLs: $($dlls.Count)"
Write-Host "[package] zip: $zip" -ForegroundColor Green
