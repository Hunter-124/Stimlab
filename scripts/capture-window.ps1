# capture-window.ps1 - launch StimLab, screenshot its window, kill it.
# DX11 swapchains don't reliably honor PrintWindow, so we foreground the window
# and BitBlt the on-screen region. computer-use can't see the uninstalled exe;
# this gives an agent a PNG to Read.
#
#   .\scripts\capture-window.ps1 -Exe build\windows\bin\StimLab.exe -Out out.png -Panel Settings
param(
    [string]$Exe = "build\windows\bin\StimLab.exe",
    [string]$Out = "stimlab.png",
    [string]$Panel = "",
    [string]$Molecule = "",
    [int]$WaitMs = 3500
)
$ErrorActionPreference = "Stop"

Add-Type @"
using System;
using System.Runtime.InteropServices;
public static class W {
    [DllImport("user32.dll")] public static extern IntPtr FindWindow(string c, string n);
    [DllImport("user32.dll")] public static extern bool SetForegroundWindow(IntPtr h);
    [DllImport("user32.dll")] public static extern bool ShowWindow(IntPtr h, int n);
    [DllImport("user32.dll")] public static extern bool SetWindowPos(IntPtr h, IntPtr after, int x, int y, int cx, int cy, uint flags);
    [DllImport("user32.dll")] public static extern bool GetClientRect(IntPtr h, out RECT r);
    [DllImport("user32.dll")] public static extern bool ClientToScreen(IntPtr h, ref POINT p);
    public struct RECT { public int L, T, R, B; }
    public struct POINT { public int X, Y; }
}
"@

if ($Panel)    { $env:STIMLAB_PANEL = $Panel }
if ($Molecule) { $env:STIMLAB_MOLECULE = $Molecule }

$proc = Start-Process -FilePath $Exe -PassThru
try {
    Start-Sleep -Milliseconds $WaitMs
    # Poll the started process for its main window handle (robust vs title match).
    $hwnd = [IntPtr]::Zero
    for ($i = 0; $i -lt 20; $i++) {
        $proc.Refresh()
        $hwnd = $proc.MainWindowHandle
        if ($hwnd -ne [IntPtr]::Zero) { break }
        Start-Sleep -Milliseconds 500
    }
    if ($hwnd -eq [IntPtr]::Zero) { throw "StimLab window not found" }
    [W]::ShowWindow($hwnd, 9) | Out-Null      # SW_RESTORE
    # Force topmost so it sits above other (non-topmost) windows - SetForegroundWindow
    # can't always steal focus, which leaves the window occluded for the BitBlt.
    $HWND_TOPMOST = [IntPtr](-1); $HWND_NOTOPMOST = [IntPtr](-2)
    $SWP = 0x40 -bor 0x1 -bor 0x2   # SHOWWINDOW | NOMOVE | NOSIZE
    [W]::SetWindowPos($hwnd, $HWND_TOPMOST, 0,0,0,0, $SWP) | Out-Null
    [W]::SetForegroundWindow($hwnd) | Out-Null
    Start-Sleep -Milliseconds 900
    $r = New-Object W+RECT
    [W]::GetClientRect($hwnd, [ref]$r) | Out-Null
    $tl = New-Object W+POINT; $tl.X = 0; $tl.Y = 0
    [W]::ClientToScreen($hwnd, [ref]$tl) | Out-Null
    $w = $r.R - $r.L; $h = $r.B - $r.T
    Add-Type -AssemblyName System.Drawing
    $bmp = New-Object System.Drawing.Bitmap($w, $h)
    $g = [System.Drawing.Graphics]::FromImage($bmp)
    $g.CopyFromScreen($tl.X, $tl.Y, 0, 0, (New-Object System.Drawing.Size($w, $h)))
    $bmp.Save((Resolve-Path -LiteralPath (Split-Path $Out -Parent)).Path + "\" + (Split-Path $Out -Leaf))
    Write-Host "[capture] wrote $Out ($w x $h)"
}
finally {
    if ($proc -and -not $proc.HasExited) { Stop-Process -Id $proc.Id -Force }
    Remove-Item Env:STIMLAB_PANEL -ErrorAction SilentlyContinue
    Remove-Item Env:STIMLAB_MOLECULE -ErrorAction SilentlyContinue
}
