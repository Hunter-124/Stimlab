# sign.ps1 - optional Authenticode signing of BioCAD.exe.
#
# BioCAD ships UNSIGNED by default. That is fine for personal use but trips
# SmartScreen ("Windows protected your PC") on first run on other machines. If you
# have a code-signing certificate, this signs the exe with signtool; WITHOUT a cert
# it is a clean no-op (exit 0) so CI stays green and the build is never blocked.
#
# Provide a certificate one of two ways (env vars):
#   $env:BIOCAD_SIGN_PFX  = path to a .pfx file   (+ $env:BIOCAD_SIGN_PFX_PASSWORD)
#   $env:BIOCAD_SIGN_SHA1 = SHA-1 thumbprint of a cert already in the Windows store
# Optional: $env:BIOCAD_SIGN_TS_URL = RFC-3161 timestamp URL (default DigiCert).
#
#   .\scripts\sign.ps1 -Exe dist\BioCAD-0.1.0-win-x64\BioCAD.exe
param(
    [Parameter(Mandatory = $true)][string]$Exe
)
$ErrorActionPreference = "Stop"
if (-not (Test-Path $Exe)) { Write-Host "[sign] exe not found: $Exe" -ForegroundColor Red; exit 1 }

$pfx  = $env:BIOCAD_SIGN_PFX
$sha1 = $env:BIOCAD_SIGN_SHA1
if (-not $pfx -and -not $sha1) {
    Write-Host "[sign] no signing certificate configured (set BIOCAD_SIGN_PFX or BIOCAD_SIGN_SHA1)." -ForegroundColor Yellow
    Write-Host "[sign] leaving BioCAD.exe UNSIGNED - SmartScreen warning on first run elsewhere is expected." -ForegroundColor Yellow
    exit 0   # signing is optional; absence of a cert is not a build failure.
}

# Locate signtool.exe (Windows SDK; not always on PATH).
$signtool = (Get-Command signtool.exe -ErrorAction SilentlyContinue).Source
if (-not $signtool) {
    $cand = Get-ChildItem "C:\Program Files (x86)\Windows Kits\10\bin\*\x64\signtool.exe" -ErrorAction SilentlyContinue |
            Sort-Object FullName -Descending | Select-Object -First 1
    if ($cand) { $signtool = $cand.FullName }
    else { throw "signtool.exe not found - install the Windows 10/11 SDK (or add it to PATH)." }
}

$ts = if ($env:BIOCAD_SIGN_TS_URL) { $env:BIOCAD_SIGN_TS_URL } else { "http://timestamp.digicert.com" }
$signArgs = @("sign", "/fd", "SHA256", "/tr", $ts, "/td", "SHA256")
if ($pfx) {
    if (-not (Test-Path $pfx)) { throw "BIOCAD_SIGN_PFX points to a missing file: $pfx" }
    $signArgs += @("/f", $pfx)
    if ($env:BIOCAD_SIGN_PFX_PASSWORD) { $signArgs += @("/p", $env:BIOCAD_SIGN_PFX_PASSWORD) }
} else {
    $signArgs += @("/sha1", $sha1)
}
$signArgs += $Exe

Write-Host "[sign] signing $Exe (signtool: $signtool)" -ForegroundColor Cyan
& $signtool @signArgs
if ($LASTEXITCODE -ne 0) { throw "signtool sign failed (exit $LASTEXITCODE)" }
& $signtool verify /pa $Exe
if ($LASTEXITCODE -ne 0) { throw "signtool verify failed (exit $LASTEXITCODE)" }
Write-Host "[sign] signed + timestamped + verified." -ForegroundColor Green
