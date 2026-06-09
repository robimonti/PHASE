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
    [string]$Source = (Join-Path $PSScriptRoot 'install-phase.ps1'),
    [string]$Output = (Join-Path $PSScriptRoot 'install-phase.exe'),
    [string]$IconFile,
    [switch]$Force
)

$ErrorActionPreference = 'Stop'

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
    description = 'PHASE Windows installer - clona repo + configura paths'
    company    = 'pyccino'
    product    = 'PHASE'
    version    = '1.1.0.0'
    noConsole  = $true
    requireAdmin = $false
    STA        = $true
}
if ($IconFile -and (Test-Path $IconFile)) {
    $ps2exeArgs.iconFile = $IconFile
}

Invoke-PS2EXE @ps2exeArgs

if (Test-Path $Output) {
    $size = (Get-Item $Output).Length / 1MB
    Write-Host "✓ Compilato: $Output ($([math]::Round($size, 2)) MB)" -ForegroundColor Green
    Write-Host ""
    Write-Host "Per distribuire come pacchetto completo (con SNAP bundled):"
    Write-Host "  1. mkdir phase-installer-package"
    Write-Host "  2. copy install-phase.exe phase-installer-package\"
    Write-Host "  3. mkdir phase-installer-package\installers"
    Write-Host "  4. copy F:\phase\installers\esa-snap_sentinel_windows-13.0.0.exe phase-installer-package\installers\"
    Write-Host "  5. Compress-Archive phase-installer-package phase-installer-v1.0.0.zip"
} else {
    throw "Compilazione fallita: $Output non creato."
}
