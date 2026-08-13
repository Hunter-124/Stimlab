# ci.ps1 - local continuous integration for BioCAD.
#
# Builds + ctests the shipping presets and produces the release zip, exactly as the
# GitHub Actions workflow (.github/workflows/ci.yml) does - so "green locally" ==
# "green in CI". Exits non-zero on the first failure (a CI gate, not best-effort).
#
#   .\scripts\ci.ps1                 # windows + windows-static: build, ctest, package
#   .\scripts\ci.ps1 -Science        # ALSO build+test windows-science-static (static curl)
#   .\scripts\ci.ps1 -Sign           # Authenticode-sign the packaged exe (needs a cert)
#   .\scripts\ci.ps1 -NoPackage      # skip the package/zip step
[CmdletBinding()]
param(
    [switch]$Science,
    [switch]$Sign,
    [switch]$NoPackage,
    [string]$Version = ""
)
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root
if (-not $Version) { $Version = & "$PSScriptRoot\version.ps1" }

# The shipping presets: the fast dynamic dev build and the static single-exe release.
# -Science adds the static build that also bundles the live agent + web tools (curl).
$presets = @("windows", "windows-static")
if ($Science) { $presets += "windows-science-static" }

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host " BioCAD CI  |  presets: $($presets -join ', ')" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan

foreach ($p in $presets) {
    Write-Host "`n[ci] >>> build + ctest: $p" -ForegroundColor Cyan
    # dev-build.ps1 enters MSVC vcvars and pins VCPKG_ROOT (the wrapper the whole repo uses).
    & "$root\scripts\dev-build.ps1" $p -Test
    if ($LASTEXITCODE -ne 0) {
        Write-Host "[ci] FAILED on preset '$p' (exit $LASTEXITCODE)" -ForegroundColor Red
        exit $LASTEXITCODE
    }
}

if (-not $NoPackage) {
    Write-Host "`n[ci] >>> packaging windows-static release zip" -ForegroundColor Cyan
    & "$root\scripts\package.ps1" -NoBuild -Preset windows-static -Version $Version
    if ($LASTEXITCODE -ne 0) { Write-Host "[ci] packaging FAILED (exit $LASTEXITCODE)" -ForegroundColor Red; exit $LASTEXITCODE }

    if ($Sign) {
        # Optional Authenticode signing. A clean no-op (exit 0) when no cert is configured.
        & "$root\scripts\sign.ps1" -Exe (Join-Path $root "dist\BioCAD-$Version-win-x64\BioCAD.exe")
        if ($LASTEXITCODE -ne 0) { Write-Host "[ci] signing FAILED (exit $LASTEXITCODE)" -ForegroundColor Red; exit $LASTEXITCODE }
    }
}

Write-Host "`n------------------------------------------------------------" -ForegroundColor Green
Write-Host " CI OK  |  presets green: $($presets -join ', ')" -ForegroundColor Green
if (-not $NoPackage) { Write-Host "   zip: $(Join-Path $root "dist\BioCAD-$Version-win-x64.zip")" }
Write-Host "------------------------------------------------------------" -ForegroundColor Green
