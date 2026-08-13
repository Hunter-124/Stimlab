# version.ps1 - single source of truth for the BioCAD version.
#
# project(BioCAD VERSION x.y.z ...) in the top-level CMakeLists.txt is authoritative;
# no script hard-codes a version number. Dot-source or invoke this to read it:
#   $ver = & "$PSScriptRoot\version.ps1"
$ErrorActionPreference = "Stop"
$root  = Split-Path -Parent $PSScriptRoot
$cmake = Join-Path $root "CMakeLists.txt"
$text  = Get-Content $cmake -Raw
if ($text -notmatch '(?m)^\s*project\(BioCAD\s*\r?\n\s*VERSION\s+([0-9]+\.[0-9]+\.[0-9]+)') {
    throw "could not parse 'project(BioCAD VERSION x.y.z' out of $cmake"
}
$Matches[1]
