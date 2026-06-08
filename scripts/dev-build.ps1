# dev-build.ps1 - reproducible StimLab build wrapper.
#
# Why this exists: vcvars64.bat OVERRIDES $env:VCPKG_ROOT to the Visual-Studio-bundled
# vcpkg, whose registry baseline differs from the user's vcpkg (and lacks some ports).
# This wrapper enters the MSVC x64 environment AND then forces VCPKG_ROOT back to the
# user's vcpkg so dependency resolution is consistent and reproducible.
#
# Usage:
#   .\scripts\dev-build.ps1                      # configure + build (windows preset)
#   .\scripts\dev-build.ps1 -Test                # configure + build + ctest
#   .\scripts\dev-build.ps1 windows-science      # Phase C preset (RDKit/curl)
#   .\scripts\dev-build.ps1 -ConfigureOnly       # configure only

param(
    [string]$Preset = "windows",
    [switch]$Test,
    [switch]$ConfigureOnly
)

$ErrorActionPreference = "Stop"

$VcpkgRoot = if ($env:STIMLAB_VCPKG_ROOT) { $env:STIMLAB_VCPKG_ROOT } else { "C:\Users\nigga\vcpkg" }
$Vcvars = "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
if (-not (Test-Path $Vcvars)) {
    $Vcvars = "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
}

$steps = "set VCPKG_ROOT=$VcpkgRoot && cmake --preset $Preset"
if (-not $ConfigureOnly) { $steps += " && cmake --build --preset $Preset" }
if ($Test) { $steps += " && ctest --preset $Preset --output-on-failure" }

Write-Host "[dev-build] preset=$Preset VCPKG_ROOT=$VcpkgRoot" -ForegroundColor Cyan
cmd /c "`"$Vcvars`" && $steps"
exit $LASTEXITCODE
