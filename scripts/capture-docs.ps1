# capture-docs.ps1 - render every panel to docs/media with the in-app --shot path.
#
# This replaces desktop screen-scraping: BioCAD.exe renders its own back buffer to a
# PNG, so it works on a CI runner with no interactive desktop session, and every
# image is exactly the same size regardless of DPI or theme metrics.
#
# Usage:
#   .\scripts\capture-docs.ps1                      # capture all panels from build\windows
#   .\scripts\capture-docs.ps1 -Preset windows-static
#   .\scripts\capture-docs.ps1 -Only Docking,Workflows
param(
    [string]$Preset  = "windows",
    [string]$OutDir  = "docs/media",
    [string[]]$Only  = @(),
    [int]$Width      = 1600,
    [int]$Height     = 1000,
    [int]$Warmup     = 120
)
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root

$exe = Join-Path $root "build\$Preset\bin\BioCAD.exe"
if (-not (Test-Path $exe)) { throw "BioCAD.exe not found at $exe - build the '$Preset' preset first." }

$outFull = Join-Path $root $OutDir
New-Item -ItemType Directory -Force -Path $outFull | Out-Null

# The panel ids are the PanelInfo ids in src/ui/AppShell.cpp; --shot-panel sets
# UiState::activePanel before the first frame. `env` forces interesting state for the
# panels that would otherwise capture empty.
$shots = @(
    @{ name = "dashboard";   panel = "Dashboard" }
    @{ name = "structure";   panel = "Structure" }
    @{ name = "input";       panel = "Input" }
    @{ name = "library";     panel = "Library" }
    @{ name = "stability";   panel = "Stability" }
    @{ name = "absorption";  panel = "Absorption" }
    @{ name = "metabolism";  panel = "Metabolism" }
    @{ name = "analog";      panel = "Analog";   env = @{ BIOCAD_ANALOG_DRAW = "1" } }
    @{ name = "compare";     panel = "Compare" }
    @{ name = "similarity";  panel = "Similarity" }
    @{ name = "docking";     panel = "Docking";  env = @{ BIOCAD_TARGET = "DAT" } }
    @{ name = "workflows";   panel = "Workflows" }
    @{ name = "legal";       panel = "Legal" }
    @{ name = "runs";        panel = "Runs" }
    @{ name = "presets";     panel = "Presets" }
    @{ name = "pkpd";        panel = "PkPd" }
    @{ name = "settings";    panel = "Settings" }
)

# Frame sequences, encoded to GIFs by the Linux CI job. The Workflows DAG and the
# docking 3D viewport are the two things a still image cannot show.
$sequences = @(
    @{ name = "workflows-dag"; panel = "Workflows"; frames = 90 }
    @{ name = "docking-3d";    panel = "Docking";   frames = 90
       env = @{ BIOCAD_TARGET = "DAT"; BIOCAD_DOCK_SCROLL3D = "1" } }
)

function Invoke-Shot {
    param($Name, $Panel, $Frames, $Env)

    $saved = @{}
    if ($Env) {
        foreach ($k in $Env.Keys) {
            $saved[$k] = [Environment]::GetEnvironmentVariable($k)
            [Environment]::SetEnvironmentVariable($k, $Env[$k])
        }
    }
    try {
        $target = Join-Path $outFull "$Name.png"
        $args = @("--shot", $target, "--shot-panel", $Panel,
                  "--shot-warmup", $Warmup, "--shot-size", $Width, $Height)
        if ($Frames -gt 1) { $args += @("--shot-frames", $Frames) }

        Write-Host "[capture] $Name ($Panel)" -ForegroundColor Cyan
        & $exe @args
        if ($LASTEXITCODE -ne 0) { throw "$Name : BioCAD.exe exited $LASTEXITCODE (3 = capture failed)" }

        # A capture that produced nothing, or a suspiciously tiny file, is a failure:
        # a 0-byte or all-black PNG in the docs is worse than no PNG at all.
        $produced = if ($Frames -gt 1) {
            Get-ChildItem $outFull -Filter "$Name-*.png"
        } else {
            Get-Item $target -ErrorAction SilentlyContinue
        }
        if (-not $produced) { throw "$Name : no PNG was written" }
        foreach ($f in @($produced)) {
            if ($f.Length -lt 4096) { throw "$($f.Name) is $($f.Length) bytes - the render did not settle" }
        }
    } finally {
        foreach ($k in $saved.Keys) { [Environment]::SetEnvironmentVariable($k, $saved[$k]) }
    }
}

$wanted = { param($n, $p) ($Only.Count -eq 0) -or ($Only -contains $n) -or ($Only -contains $p) }

foreach ($s in $shots) {
    if (-not (& $wanted $s.name $s.panel)) { continue }
    Invoke-Shot -Name $s.name -Panel $s.panel -Frames 1 -Env $s.env
}
foreach ($s in $sequences) {
    if (-not (& $wanted $s.name $s.panel)) { continue }
    Invoke-Shot -Name $s.name -Panel $s.panel -Frames $s.frames -Env $s.env
}

$count = (Get-ChildItem $outFull -Filter *.png).Count
Write-Host "[capture] $count PNG(s) in $OutDir" -ForegroundColor Green
