# sign.ps1 - optional Authenticode signing of StimLab.exe.
#
# StimLab ships UNSIGNED by default. That is fine for personal use but trips
# SmartScreen ("Windows protected your PC") on first run on other machines. If you
# have a code-signing certificate, this signs the exe with signtool; WITHOUT a cert
# it is a clean no-op (exit 0) so CI stays green and the build is never blocked.
#
# Provide a certificate one of two ways (env vars):
#   $env:STIMLAB_SIGN_PFX  = path to a .pfx file   (+ $env:STIMLAB_SIGN_PFX_PASSWORD)
#   $env:STIMLAB_SIGN_SHA1 = SHA-1 thumbprint of a cert already in the Windows store
# Optional: $env:STIMLAB_SIGN_TS_URL = RFC-3161 timestamp URL (default DigiCert).
#
#   .\scripts\sign.ps1 -Exe dist\StimLab-0.1.0-win-x64\StimLab.exe
param(
    [Parameter(Mandatory = $true)][string]$Exe
)
$ErrorActionPreference = "Stop"
if (-not (Test-Path $Exe)) { Write-Host "[sign] exe not found: $Exe" -ForegroundColor Red; exit 1 }

$pfx  = $env:STIMLAB_SIGN_PFX
$sha1 = $env:STIMLAB_SIGN_SHA1
if (-not $pfx -and -not $sha1) {
    Write-Host "[sign] no signing certificate configured (set STIMLAB_SIGN_PFX or STIMLAB_SIGN_SHA1)." -ForegroundColor Yellow
    Write-Host "[sign] leaving StimLab.exe UNSIGNED - SmartScreen warning on first run elsewhere is expected." -ForegroundColor Yellow
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

$ts = if ($env:STIMLAB_SIGN_TS_URL) { $env:STIMLAB_SIGN_TS_URL } else { "http://timestamp.digicert.com" }
$signArgs = @("sign", "/fd", "SHA256", "/tr", $ts, "/td", "SHA256")
if ($pfx) {
    if (-not (Test-Path $pfx)) { throw "STIMLAB_SIGN_PFX points to a missing file: $pfx" }
    $signArgs += @("/f", $pfx)
    if ($env:STIMLAB_SIGN_PFX_PASSWORD) { $signArgs += @("/p", $env:STIMLAB_SIGN_PFX_PASSWORD) }
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
