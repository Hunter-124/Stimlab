# build.ps1 - one-command build for BioCAD.
#
# Wraps the MSVC/vcpkg build so you don't have to remember presets or vcvars. Two modes:
#
#   .\build.ps1                 # fast DEV build (dynamic) + run the test suite
#   .\build.ps1 -Release        # static single-exe RELEASE build + tests + packaged zip
#   .\build.ps1 -Release -Clean # ... wiping the build dir first (from-scratch)
#   .\build.ps1 -NoTest         # skip ctest
#   .\build.ps1 -Science        # dev build WITH the live LLM + web tools (libcurl)
#   .\build.ps1 -Release -Science # static single-exe release that ALSO bundles the
#                                 # live agent + web tools (curl[ssl], Schannel, static)
#
# -Release produces dist\BioCAD-<Version>-win-x64\BioCAD.exe (a self-contained exe
# needing no DLLs / no VC++ redist) and a matching .zip. Adding -Science to -Release
# makes that single exe also carry the live Anthropic provider + web_search/web_fetch.
[CmdletBinding()]
param(
    [switch]$Release,
    [switch]$Science,
    [switch]$NoTest,
    [switch]$Clean,
    [string]$Version = ""
)
$ErrorActionPreference = "Stop"
$root = $PSScriptRoot
Set-Location $root
if (-not $Version) { $Version = & "$root\scripts\version.ps1" }

# Pick the preset. -Release wins (static, packaged); -Release -Science is the static
# release that also bundles the live agent + web tools; else dev (dynamic), optionally science.
$preset = if ($Release -and $Science) { "windows-science-static" }
          elseif ($Release)           { "windows-static" }
          elseif ($Science)           { "windows-science" }
          else                        { "windows" }

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host " BioCAD build  |  preset=$preset  test=$(-not $NoTest)  release=$Release" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

if ($Clean) {
    $bd = Join-Path $root "build\$preset"
    if (Test-Path $bd) {
        Write-Host "[build] cleaning $bd ..." -ForegroundColor DarkYellow
        Remove-Item -Recurse -Force $bd
    }
}

# Configure + build (+ ctest unless -NoTest). dev-build.ps1 enters vcvars and pins VCPKG_ROOT.
$devArgs = @($preset)
if (-not $NoTest) { $devArgs += "-Test" }
& "$root\scripts\dev-build.ps1" @devArgs
if ($LASTEXITCODE -ne 0) { Write-Host "[build] FAILED (exit $LASTEXITCODE)" -ForegroundColor Red; exit $LASTEXITCODE }

$exe = Join-Path $root "build\$preset\bin\BioCAD.exe"
if (-not (Test-Path $exe)) { Write-Host "[build] BioCAD.exe not produced" -ForegroundColor Red; exit 1 }
$exeMB = "{0:N1} MB" -f ((Get-Item $exe).Length / 1MB)

if ($Release) {
    Write-Host "[build] packaging release ..." -ForegroundColor Cyan
    & "$root\scripts\package.ps1" -NoBuild -Preset $preset -Version $Version
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}

Write-Host "------------------------------------------------------------" -ForegroundColor Green
Write-Host " BUILD OK" -ForegroundColor Green
Write-Host "   exe : $exe ($exeMB)"
if ($Release) {
    Write-Host "   dist: $(Join-Path $root "dist\BioCAD-$Version-win-x64")"
    Write-Host "   zip : $(Join-Path $root "dist\BioCAD-$Version-win-x64.zip")"
}
Write-Host "------------------------------------------------------------" -ForegroundColor Green
