# Compila install-phase.ps1 in install-phase.exe via PS2EXE.
#
# Una tantum, installa PS2EXE:
#   Install-Module -Name ps2exe -Scope CurrentUser -Force
#
# Poi:
#   powershell -ExecutionPolicy Bypass -File compile-to-exe.ps1
#
# Output: install-phase.exe accanto a questo script.
#
# Nota: il .exe NON include l'installer SNAP bundled (~500 MB).
# Per distribuirlo come pacchetto:
#   1. Compila .exe con questo script.
#   2. Crea un .zip con dentro:
#        install-phase.exe
#        installers\esa-snap_sentinel_windows-13.0.0.exe   (copialo da F:\phase\installers\)
#   3. L'utente finale estrae lo zip e fa doppio click su install-phase.exe.
#      Lo script cerca l'installer SNAP in .\installers\ accanto a sé.

[CmdletBinding()]
param(
    [string]$Source,
    [string]$Output,
    [string]$IconFile,
    [string]$DefaultPhaseBranch,
    [switch]$Force
)

$ErrorActionPreference = 'Stop'

# Windows PowerShell can evaluate parameter default expressions before
# $PSScriptRoot has been populated. Resolve defaults only after param() so the
# script works consistently from -File, relative paths and older PS5 hosts.
$scriptPath = $MyInvocation.MyCommand.Path
$scriptDir = if (-not [string]::IsNullOrWhiteSpace($PSScriptRoot)) {
    $PSScriptRoot
} elseif (-not [string]::IsNullOrWhiteSpace($scriptPath)) {
    Split-Path -Parent $scriptPath
} else {
    (Get-Location).Path
}
if ([string]::IsNullOrWhiteSpace($Source)) {
    $Source = Join-Path $scriptDir 'install-phase.ps1'
}
if ([string]::IsNullOrWhiteSpace($Output)) {
    $Output = Join-Path $scriptDir 'install-phase.exe'
}
if ([string]::IsNullOrWhiteSpace($IconFile)) {
    $IconFile = Join-Path $scriptDir 'PHASE.ico'
}

if (-not (Test-Path $Source)) {
    throw "Source non trovato: $Source"
}

if (-not (Get-Module -ListAvailable -Name ps2exe)) {
    Write-Host "Modulo ps2exe non installato. Lo installo nell'utente corrente..."
    Install-Module -Name ps2exe -Scope CurrentUser -Force -AllowClobber
}

Import-Module ps2exe

if ((Test-Path $Output) -and -not $Force) {
    $ans = Read-Host "$Output esiste già. Sovrascrivo? (s/N)"
    if ($ans -notin @('s','S','y','Y')) {
        Write-Host "Compilazione annullata."
        exit 0
    }
}

Write-Host "Compilazione $Source -> $Output ..." -ForegroundColor Cyan

$ps2exeArgs = @{
    inputFile  = $Source
    outputFile = $Output
    title      = 'PHASE Installer'
    description = 'PHASE 7 unified hub installer'
    company    = 'Roberto Monti and pyccino'
    product    = 'PHASE'
    version    = '7.0.0.0'
    noConsole  = $true
    requireAdmin = $false
    STA        = $true
}
if ($DefaultPhaseBranch -match '^v(7\.\d+\.\d+)$') {
    $ps2exeArgs.version = "$($Matches[1]).0"
}
if ($IconFile -and (Test-Path $IconFile)) {
    $ps2exeArgs.iconFile = $IconFile
}

$temporarySource = $null
try {
    $scriptText = Get-Content -LiteralPath $Source -Raw
    if (-not [string]::IsNullOrWhiteSpace($DefaultPhaseBranch)) {
        if ($DefaultPhaseBranch -notmatch '^[A-Za-z0-9_./-]+$') {
            throw "Branch name is not safe for embedding: $DefaultPhaseBranch"
        }
        $marker = "[string]`$PhaseBranch = 'main'"
        if (-not $scriptText.Contains($marker)) {
            throw "Expected default branch declaration not found in $Source"
        }
        $scriptText = $scriptText.Replace($marker, "[string]`$PhaseBranch = '$DefaultPhaseBranch'")
        if ($DefaultPhaseBranch -match '^v7\.\d+\.\d+$') {
            $scriptText = $scriptText.Replace('v7.0.0 preview', $DefaultPhaseBranch)
        }
        Write-Host "Embedded PHASE branch: $DefaultPhaseBranch"
    }
    $logoPath = Join-Path $scriptDir 'PHASE_logo.png'
    if (-not (Test-Path -LiteralPath $logoPath)) {
        throw "Installer logo not found: $logoPath"
    }
    $logoMarker = "`$Script:EmbeddedLogoBase64 = ''"
    if (-not $scriptText.Contains($logoMarker)) {
        throw "Expected embedded logo marker not found in $Source"
    }
    $logoBase64 = [Convert]::ToBase64String([IO.File]::ReadAllBytes($logoPath))
    $scriptText = $scriptText.Replace($logoMarker, "`$Script:EmbeddedLogoBase64 = '$logoBase64'")
    $temporarySource = Join-Path $env:TEMP ("phase-installer-" + [guid]::NewGuid().ToString('N') + '.ps1')
    Set-Content -LiteralPath $temporarySource -Value $scriptText -Encoding UTF8
    $ps2exeArgs.inputFile = $temporarySource
    Invoke-PS2EXE @ps2exeArgs
} finally {
    if ($temporarySource -and (Test-Path -LiteralPath $temporarySource)) {
        Remove-Item -LiteralPath $temporarySource -Force
    }
}

if (Test-Path $Output) {
    $size = (Get-Item $Output).Length / 1MB
    Write-Host "✓ Compilato: $Output ($([math]::Round($size, 2)) MB)" -ForegroundColor Green
    Write-Host ""
    Write-Host "Per distribuire come pacchetto completo (con SNAP bundled):"
    Write-Host "  1. mkdir phase-installer-package"
    Write-Host "  2. copy install-phase.exe phase-installer-package\"
    Write-Host "  3. mkdir phase-installer-package\installers"
    Write-Host "  4. copy F:\phase\installers\esa-snap_sentinel_windows-13.0.0.exe phase-installer-package\installers\"
    Write-Host "  5. Compress-Archive phase-installer-package phase-installer-v7.zip"
} else {
    throw "Compilazione fallita: $Output non creato."
}
