# make-icon.ps1 - generate the StimLab app icon (src/app/StimLab.ico).
#
# Draws the brand mark (a molecule ring on an indigo->violet gradient squircle) at
# several sizes with GDI+, then packs them into a multi-resolution .ico (32-bit BGRA
# DIB entries + AND mask, the format LoadIcon/Explorer expect). Re-run to regenerate.
#
#   .\scripts\make-icon.ps1
param(
    [string]$Out = "$PSScriptRoot\..\src\app\StimLab.ico"
)
$ErrorActionPreference = "Stop"
Add-Type -AssemblyName System.Drawing

# Brand palette (matches src/ui/Theme.h).
$indigo = [System.Drawing.Color]::FromArgb(255, 110, 98, 246)
$violet = [System.Drawing.Color]::FromArgb(255, 192, 132, 252)
$white  = [System.Drawing.Color]::FromArgb(255, 245, 247, 255)

function New-RoundedPath([int]$s, [int]$r) {
    $p = New-Object System.Drawing.Drawing2D.GraphicsPath
    $d = $r * 2
    $p.AddArc(0, 0, $d, $d, 180, 90)
    $p.AddArc($s - $d, 0, $d, $d, 270, 90)
    $p.AddArc($s - $d, $s - $d, $d, $d, 0, 90)
    $p.AddArc(0, $s - $d, $d, $d, 90, 90)
    $p.CloseFigure()
    return $p
}

function New-LogoBitmap([int]$s) {
    $bmp = New-Object System.Drawing.Bitmap($s, $s, [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
    $g = [System.Drawing.Graphics]::FromImage($bmp)
    $g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
    $g.PixelOffsetMode = [System.Drawing.Drawing2D.PixelOffsetMode]::HighQuality
    $g.Clear([System.Drawing.Color]::Transparent)

    # Gradient squircle background.
    $rect = New-Object System.Drawing.Rectangle(0, 0, $s, $s)
    $radius = [int]($s * 0.22)
    $path = New-RoundedPath $s $radius
    $brush = New-Object System.Drawing.Drawing2D.LinearGradientBrush($rect, $indigo, $violet, 50.0)
    $g.FillPath($brush, $path)

    # Molecule ring (white) centred a touch low so an optional tail has room.
    $cx = $s * 0.50
    $cy = $s * 0.53
    $r  = $s * 0.255
    $pw = [Math]::Max(1.6, $s * 0.072)
    $pen = New-Object System.Drawing.Pen($white, $pw)
    $pen.LineJoin = [System.Drawing.Drawing2D.LineJoin]::Round
    $pen.StartCap = [System.Drawing.Drawing2D.LineCap]::Round
    $pen.EndCap   = [System.Drawing.Drawing2D.LineCap]::Round
    $wb = New-Object System.Drawing.SolidBrush($white)

    $verts = @()
    for ($i = 0; $i -lt 6; $i++) {
        $a = [Math]::PI / 180.0 * (-90 + $i * 60)
        $verts += , (New-Object System.Drawing.PointF([single]($cx + $r * [Math]::Cos($a)), [single]($cy + $r * [Math]::Sin($a))))
    }

    # Optional substituent tail off the upper-right vertex (skip at tiny sizes).
    if ($s -ge 32) {
        $v = $verts[1]
        $dirx = ($v.X - $cx); $diry = ($v.Y - $cy)
        $len = [Math]::Sqrt($dirx * $dirx + $diry * $diry)
        $ux = $dirx / $len; $uy = $diry / $len
        $tip = New-Object System.Drawing.PointF([single]($v.X + $ux * $r * 0.62), [single]($v.Y + $uy * $r * 0.62))
        $g.DrawLine($pen, $v, $tip)
        $dr2 = $s * 0.062
        $g.FillEllipse($wb, [single]($tip.X - $dr2), [single]($tip.Y - $dr2), [single]($dr2 * 2), [single]($dr2 * 2))
    }

    # Hexagon outline + inner aromatic circle + vertex dots.
    $g.DrawPolygon($pen, [System.Drawing.PointF[]]$verts)
    $ir = $r * 0.46
    $ipw = [Math]::Max(1.0, $s * 0.034)
    $ipen = New-Object System.Drawing.Pen($white, $ipw)
    $g.DrawEllipse($ipen, [single]($cx - $ir), [single]($cy - $ir), [single]($ir * 2), [single]($ir * 2))
    $dr = [Math]::Max(1.3, $s * 0.058)
    foreach ($v in $verts) {
        $g.FillEllipse($wb, [single]($v.X - $dr), [single]($v.Y - $dr), [single]($dr * 2), [single]($dr * 2))
    }

    $g.Dispose()
    return $bmp
}

# Build a 32bpp BGRA DIB (BITMAPINFOHEADER + XOR rows + AND mask), bottom-up.
function Get-Dib([System.Drawing.Bitmap]$bmp) {
    $w = $bmp.Width; $h = $bmp.Height
    $ms = New-Object System.IO.MemoryStream
    $bw = New-Object System.IO.BinaryWriter -ArgumentList $ms
    $bw.Write([int]40); $bw.Write([int]$w); $bw.Write([int]($h * 2))
    $bw.Write([uint16]1); $bw.Write([uint16]32)
    $bw.Write([int]0); $bw.Write([int]0); $bw.Write([int]0); $bw.Write([int]0); $bw.Write([int]0); $bw.Write([int]0)
    for ($y = $h - 1; $y -ge 0; $y--) {
        for ($x = 0; $x -lt $w; $x++) {
            $c = $bmp.GetPixel($x, $y)
            $bw.Write([byte]$c.B); $bw.Write([byte]$c.G); $bw.Write([byte]$c.R); $bw.Write([byte]$c.A)
        }
    }
    $rowBytes = [int]([Math]::Floor(($w + 31) / 32) * 4)
    $zero = New-Object byte[] $rowBytes
    for ($y = 0; $y -lt $h; $y++) { $bw.Write($zero, 0, $rowBytes) }
    $bw.Flush()
    return $ms.ToArray()
}

$sizes = @(16, 24, 32, 48, 64, 128, 256)
$imgs = @()
foreach ($s in $sizes) {
    $bmp = New-LogoBitmap $s
    $imgs += , @{ w = $s; h = $s; bytes = (Get-Dib $bmp) }
    $bmp.Dispose()
}

$icoMs = New-Object System.IO.MemoryStream
$icoBw = New-Object System.IO.BinaryWriter -ArgumentList $icoMs
$icoBw.Write([uint16]0); $icoBw.Write([uint16]1); $icoBw.Write([uint16]$imgs.Count)
$offset = 6 + 16 * $imgs.Count
foreach ($e in $imgs) {
    $icoBw.Write([byte]($(if ($e.w -ge 256) { 0 } else { $e.w })))
    $icoBw.Write([byte]($(if ($e.h -ge 256) { 0 } else { $e.h })))
    $icoBw.Write([byte]0); $icoBw.Write([byte]0)
    $icoBw.Write([uint16]1); $icoBw.Write([uint16]32)
    $icoBw.Write([int]$e.bytes.Length); $icoBw.Write([int]$offset)
    $offset += $e.bytes.Length
}
foreach ($e in $imgs) { $icoBw.Write($e.bytes, 0, $e.bytes.Length) }
$icoBw.Flush()

$resolved = [System.IO.Path]::GetFullPath($Out)
[System.IO.File]::WriteAllBytes($resolved, $icoMs.ToArray())
Write-Host "[make-icon] wrote $resolved ($([Math]::Round((Get-Item $resolved).Length/1KB,1)) KB, $($imgs.Count) sizes)"
