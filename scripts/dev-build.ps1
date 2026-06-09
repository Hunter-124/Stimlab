# dev-build.ps1 - reproducible StimLab build wrapper.
#
# Auto-detects the toolchain so it works across machines (no hardcoded user paths):
#   * Visual Studio (vcvars64.bat) via vswhere, with sensible fallbacks.
#   * A bootstrapped vcpkg via $STIMLAB_VCPKG_ROOT / $VCPKG_ROOT / common locations.
#   * The VS-bundled CMake + Ninja are put on PATH (the Ninja generator needs ninja,
#     and plain vcvars64.bat does not add the CommonExtensions CMake/Ninja).
#
# Why VCPKG_ROOT is forced AFTER vcvars: vcvars64.bat overrides $env:VCPKG_ROOT to the
# Visual-Studio-bundled vcpkg, whose registry baseline differs (and lacks some ports).
# This wrapper enters the MSVC x64 environment AND then pins VCPKG_ROOT back to the
# standalone vcpkg so dependency resolution stays consistent and reproducible.
#
# Usage:
#   .\scripts\dev-build.ps1                      # configure + build (windows preset)
#   .\scripts\dev-build.ps1 -Test                # configure + build + ctest
#   .\scripts\dev-build.ps1 windows-science      # Phase C preset (RDKit/curl)
#   .\scripts\dev-build.ps1 -ConfigureOnly       # configure only
#
# Overrides (env): STIMLAB_VCVARS (path to vcvars64.bat), STIMLAB_VCPKG_ROOT (vcpkg dir).

param(
    [string]$Preset = "windows",
    [switch]$Test,
    [switch]$ConfigureOnly
)

$ErrorActionPreference = "Stop"

function Find-Vcvars {
    if ($env:STIMLAB_VCVARS -and (Test-Path $env:STIMLAB_VCVARS)) { return $env:STIMLAB_VCVARS }
    $vswhere = "C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe"
    if (Test-Path $vswhere) {
        $vs = & $vswhere -latest -products * `
            -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 `
            -property installationPath
        if ($vs) {
            $vc = Join-Path $vs "VC\Auxiliary\Build\vcvars64.bat"
            if (Test-Path $vc) { return $vc }
        }
    }
    $fallbacks = @(
        "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat",
        "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat",
        "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvars64.bat",
        "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvars64.bat",
        "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
    )
    foreach ($f in $fallbacks) { if (Test-Path $f) { return $f } }
    throw "Could not find vcvars64.bat. Install the Visual Studio 'Desktop development with C++' workload, or set STIMLAB_VCVARS."
}

function Find-VcpkgRoot {
    $candidates = @()
    if ($env:STIMLAB_VCPKG_ROOT) { $candidates += $env:STIMLAB_VCPKG_ROOT }
    if ($env:VCPKG_ROOT)         { $candidates += $env:VCPKG_ROOT }
    $candidates += @(
        (Join-Path $env:USERPROFILE "vcpkg"),
        "C:\vcpkg",
        "C:\dev\vcpkg",
        "C:\Users\Yakub\Desktop\NEWGEORGE\vcpkg"
    )
    foreach ($c in $candidates) {
        if ($c -and (Test-Path (Join-Path $c "scripts\buildsystems\vcpkg.cmake"))) {
            return (Resolve-Path $c).Path
        }
    }
    throw "Could not find a bootstrapped vcpkg (checked STIMLAB_VCPKG_ROOT, VCPKG_ROOT and common paths). Clone + bootstrap vcpkg, then set STIMLAB_VCPKG_ROOT."
}

$Vcvars    = Find-Vcvars
$VcpkgRoot = Find-VcpkgRoot

# Tool dirs to put on PATH BEFORE vcvars runs:
#   * The VS Installer dir holds vswhere.exe, which vcvars probes by bare name - without
#     it on PATH, VS 2026's vcvars prints a harmless "'vswhere.exe' is not recognized".
#   * The VS-bundled CMake + Ninja (under CommonExtensions) are needed by the Ninja
#     generator and are NOT added by vcvars itself. Locate them off the VS root (3 levels
#     up from the vcvars dir).
$VsRoot       = (Get-Item $Vcvars).Directory.Parent.Parent.Parent.FullName
$InstallerDir = "C:\Program Files (x86)\Microsoft Visual Studio\Installer"
$CMakeBin     = Join-Path $VsRoot "Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin"
$NinjaBin     = Join-Path $VsRoot "Common7\IDE\CommonExtensions\Microsoft\CMake\Ninja"
$prePath = @($InstallerDir, $NinjaBin, $CMakeBin) | Where-Object { Test-Path $_ }
$prefix  = ($prePath -join ";")

# Quote each cmd `set "NAME=value"` so a trailing space before && is NOT captured into
# the value (an unquoted `set VAR=x && cmd` stores "x " - a trailing space that corrupts
# the path and breaks a from-scratch configure; a populated cache hides it).
# PATH is set first (vcvars then augments it); VCPKG_ROOT is pinned AFTER vcvars because
# vcvars overrides it to the VS-bundled vcpkg.
$steps = ""
if ($prefix) { $steps += "set `"PATH=$prefix;%PATH%`" && " }
$steps += "`"$Vcvars`" && set `"VCPKG_ROOT=$VcpkgRoot`" && cmake --preset $Preset"
if (-not $ConfigureOnly) { $steps += " && cmake --build --preset $Preset" }
if ($Test) { $steps += " && ctest --preset $Preset --output-on-failure" }

Write-Host "[dev-build] preset      = $Preset" -ForegroundColor Cyan
Write-Host "[dev-build] vcvars      = $Vcvars" -ForegroundColor DarkGray
Write-Host "[dev-build] VCPKG_ROOT  = $VcpkgRoot" -ForegroundColor DarkGray

# Native tools (vcvars/cmake/cl) write progress and warnings to stderr; do not let
# PowerShell promote those to a terminating error (which would mask the real exit code,
# especially when this script's output is being captured to a log).
$ErrorActionPreference = "Continue"
cmd /c $steps
exit $LASTEXITCODE
