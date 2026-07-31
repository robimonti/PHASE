param(
    [string]$PhaseRoot = $PSScriptRoot
)

$ErrorActionPreference = 'Stop'
$PhaseRoot = [System.IO.Path]::GetFullPath($PhaseRoot)
$StampsRoot = Join-Path $PhaseRoot 'StaMPS'
$TrainRoot = Join-Path $PhaseRoot 'TRAIN'
$BinariesUrl = 'https://github.com/pyccino/StaMPS/releases/download/windows-port-bins-v1/stamps-win64-binaries.zip'

function Write-Step {
    param([string]$Message)
    Write-Host "[PHASE beta] $Message" -ForegroundColor Cyan
}

function Invoke-Git {
    param([string[]]$Arguments)
    & $script:GitExe @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "git failed (exit code $LASTEXITCODE): $($Arguments -join ' ')"
    }
}

function Sync-Repository {
    param(
        [string]$Repository,
        [string]$Branch,
        [string]$Destination
    )
    if (Test-Path (Join-Path $Destination '.git')) {
        Write-Step "Updating $(Split-Path $Destination -Leaf)..."
        Invoke-Git @('-C',$Destination,'fetch','origin',$Branch)
        Invoke-Git @('-C',$Destination,'checkout',$Branch)
        Invoke-Git @('-C',$Destination,'pull','--ff-only','origin',$Branch)
        return
    }
    if (Test-Path $Destination) {
        $contents = @(Get-ChildItem -LiteralPath $Destination -Force)
        if ($contents.Count -gt 0) {
            throw "$Destination exists but is not a git clone. Rename/remove that dependency folder and retry."
        }
    }
    Write-Step "Cloning $(Split-Path $Destination -Leaf)..."
    Invoke-Git @('clone','--branch',$Branch,'--single-branch',$Repository,$Destination)
}

function Install-StampsBinaries {
    param([string]$Destination)

    $routes = @{
        'calamp.exe'       = (Join-Path $Destination 'bin')
        'cpxsum.exe'       = (Join-Path $Destination 'bin')
        'pscphase.exe'     = (Join-Path $Destination 'bin')
        'pscdem.exe'       = (Join-Path $Destination 'bin')
        'psclonlat.exe'    = (Join-Path $Destination 'bin')
        'selpsc_patch.exe' = (Join-Path $Destination 'bin')
        'selsbc_patch.exe' = (Join-Path $Destination 'bin')
        'triangle.exe'     = (Join-Path $Destination 'external\triangle\bin')
        'snaphu.exe'       = (Join-Path $Destination 'external\snaphu\bin')
    }
    $missing = @($routes.Keys | Where-Object {
        -not (Test-Path (Join-Path $routes[$_] $_))
    })
    if ($missing.Count -eq 0) {
        Write-Step 'All 9 native StaMPS executables are already present.'
        return
    }

    Write-Step "Downloading the StaMPS Windows runtime ($($missing.Count) executable(s) missing)..."
    $temporaryZip = Join-Path $env:TEMP "phase-stamps-binaries-$([guid]::NewGuid().ToString('N')).zip"
    try {
        Invoke-WebRequest -Uri $BinariesUrl -OutFile $temporaryZip -UseBasicParsing
        Add-Type -AssemblyName System.IO.Compression.FileSystem
        $archive = [System.IO.Compression.ZipFile]::OpenRead($temporaryZip)
        try {
            foreach ($entry in $archive.Entries) {
                if (-not $routes.ContainsKey($entry.Name)) { continue }
                $destinationFolder = $routes[$entry.Name]
                New-Item -ItemType Directory -Path $destinationFolder -Force | Out-Null
                $destinationFile = Join-Path $destinationFolder $entry.Name
                if (Test-Path $destinationFile) {
                    Remove-Item -LiteralPath $destinationFile -Force
                }
                [System.IO.Compression.ZipFileExtensions]::ExtractToFile(
                    $entry, $destinationFile
                )
                Unblock-File -LiteralPath $destinationFile -ErrorAction SilentlyContinue
            }
        } finally {
            $archive.Dispose()
        }
    } finally {
        Remove-Item -LiteralPath $temporaryZip -Force -ErrorAction SilentlyContinue
    }

    $stillMissing = @($routes.Keys | Where-Object {
        -not (Test-Path (Join-Path $routes[$_] $_))
    })
    if ($stillMissing.Count -gt 0) {
        throw "StaMPS runtime is incomplete after extraction: $($stillMissing -join ', ')"
    }
    Write-Step 'All 9 native StaMPS executables are ready.'
}

if (-not (Test-Path (Join-Path $PhaseRoot 'PHASE_Preprocessing\PHASE_StaMPS_beta.m'))) {
    throw "PHASE standalone beta was not found in $PhaseRoot. Run this script from the extracted package root."
}

$gitCommand = Get-Command git.exe -ErrorAction SilentlyContinue
if (-not $gitCommand) {
    $gitCommand = Get-Command git -ErrorAction SilentlyContinue
}
if (-not $gitCommand) {
    throw 'Git is required for the dependency bootstrap. Install Git for Windows, reopen PowerShell and retry.'
}
$script:GitExe = $gitCommand.Source

Write-Step "Preparing external runtimes in $PhaseRoot"
Sync-Repository -Repository 'https://github.com/pyccino/StaMPS.git' `
    -Branch 'master' -Destination $StampsRoot
Sync-Repository -Repository 'https://github.com/pyccino/TRAIN.git' `
    -Branch 'main' -Destination $TrainRoot
Install-StampsBinaries -Destination $StampsRoot

$required = @(
    (Join-Path $StampsRoot 'matlab\stamps.m'),
    (Join-Path $StampsRoot 'matlab\setparm.m'),
    (Join-Path $StampsRoot 'external\snaphu\bin\snaphu.exe')
)
$missingRequired = @($required | Where-Object { -not (Test-Path $_) })
if ($missingRequired.Count -gt 0) {
    throw "Dependency verification failed: $($missingRequired -join ', ')"
}

Write-Host ''
Write-Host 'PHASE beta dependencies are ready.' -ForegroundColor Green
Write-Host "StaMPS: $StampsRoot"
Write-Host "TRAIN:  $TrainRoot"
Write-Host ''
Write-Host 'Next: restart MATLAB, cd to this folder, run addpath(genpath(pwd)) and launch PHASE_Preprocessing_beta.'
