# =============================================================================
# PHASE Windows Installer - wizard end-to-end per il porting Windows di PHASE.
#
# Cosa fa, nell'ordine:
#   1. Welcome (logo + intro).
#   2. Verifica MATLAB (auto-detect + override).
#   3. Verifica SNAP   (auto-detect + override + lancio installer bundled).
#   4. Verifica/installa Python 3.11+ (silent install da python.org se assente).
#   5. Sceglie cartella destinazione (default: Desktop\PHASE).
#   6. Clona PHASE, StaMPS, TRAIN sotto <dest>\engine\ (motore nascosto).
#   7. Scarica i binari nativi StaMPS precompilati (stamps-win64-binaries.zip).
#   8. Configura tutto: MATLAB_EXE, %APPDATA%\PHASE\python.txt, savepath MATLAB.
#   9. Crea nella root <dest> i collegamenti .lnk ai 3 .mlapp + un README,
#      così l'utente avvia l'app senza entrare in engine\.
#
# Usage (sorgente):
#   powershell -ExecutionPolicy Bypass -File install-phase.ps1
#
# Per compilare in .exe distribuibile vedere compile-to-exe.ps1.
# =============================================================================

[CmdletBinding()]
param(
    [string]$DefaultInstallDir = "$env:USERPROFILE\Desktop",
    [switch]$DryRun
)

$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName PresentationFramework
Add-Type -AssemblyName System.Windows.Forms

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------
$Script:PhaseRepo  = 'https://github.com/pyccino/PHASE.git'
$Script:PhaseBranch = 'main'   # rinominato da windows-port/main -> main (commit dffa675)
$Script:StampsRepo = 'https://github.com/pyccino/StaMPS.git'   # fork con TS picker + GUI fixes
$Script:StampsBranch = 'master'
$Script:TrainRepo  = 'https://github.com/pyccino/TRAIN.git'
$Script:TrainBranch = 'main'
$Script:PythonUrl  = 'https://www.python.org/ftp/python/3.11.9/python-3.11.9-amd64.exe'
$Script:PythonMinMajor = 3
$Script:PythonMinMinor = 11
$Script:GitPortableUrl = 'https://github.com/git-for-windows/git/releases/download/v2.45.2.windows.1/PortableGit-2.45.2-64-bit.7z.exe'
$Script:GmtZipUrl = 'https://github.com/GenericMappingTools/gmt/releases/download/6.6.0/gmt-6.6.0-win64.zip'
$Script:SnapInstallerUrl = 'https://download.esa.int/step/snap/13.0/installers/esa-snap_sentinel_windows-13.0.0.exe'

# Resolve script directory for finding bundled SNAP installer.
# Tre casi:
#   1. Sorgente .ps1 normale          -> $PSScriptRoot
#   2. Sorgente . sourced             -> $MyInvocation.MyCommand.Path
#   3. Compilato in .exe con PS2EXE   -> entrambi sopra sono null, fallback
#      al path del processo corrente.
$Script:ScriptDir = if ($PSScriptRoot) {
    $PSScriptRoot
} elseif ($MyInvocation.MyCommand.Path) {
    Split-Path -Parent $MyInvocation.MyCommand.Path
} else {
    Split-Path -Parent ([System.Diagnostics.Process]::GetCurrentProcess().MainModule.FileName)
}
$Script:BundledSnapPath = Join-Path $Script:ScriptDir 'installers\esa-snap_sentinel_windows-13.0.0.exe'

# State accumulated across wizard pages
$Script:State = @{
    MatlabExe   = $null
    SnapGpt     = $null
    PythonExe   = $null
    PythonVersion = $null
    InstallDir  = $DefaultInstallDir
    GitExe      = $null
}

# -----------------------------------------------------------------------------
# Detection helpers
# -----------------------------------------------------------------------------

# MATLAB: glob Program Files, registry HKLM Mathworks, PATH lookup.
# Returns absolute path to matlab.exe or $null.
function Find-Matlab {
    if ($env:MATLAB_EXE -and (Test-Path $env:MATLAB_EXE)) {
        return $env:MATLAB_EXE
    }

    $glob = Get-ChildItem -Path 'C:\Program Files\MATLAB\R*\bin\matlab.exe' -ErrorAction SilentlyContinue |
        Sort-Object FullName -Descending | Select-Object -First 1
    if ($glob) { return $glob.FullName }

    $reg = Get-ItemProperty -Path 'HKLM:\SOFTWARE\Mathworks\MATLAB\*' -ErrorAction SilentlyContinue |
        ForEach-Object {
            if ($_.MATLABROOT) { Join-Path $_.MATLABROOT 'bin\matlab.exe' }
        } |
        Where-Object { $_ -and (Test-Path $_) } |
        Sort-Object -Descending | Select-Object -First 1
    if ($reg) { return $reg }

    $cmd = Get-Command matlab.exe -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }

    return $null
}

# Verifica quali toolbox MATLAB sono installate ispezionando le sub-directory
# di <MATLABROOT>\toolbox\. Molto piu' veloce di matlab -batch "v=ver" (<100ms
# vs 30s di startup MATLAB). Affidabile: ogni MathWorks toolbox crea sempre
# la sua sub-folder dedicata con quel nome canonico.
#
# Returns hashtable @{ 'Mapping Toolbox' = $true/$false; ... }.
function Find-MatlabToolboxes {
    param([Parameter(Mandatory)] [string]$MatlabExe)

    # <MATLABROOT> = parent di bin\ = parent di matlab.exe parent
    $matlabRoot = Split-Path -Parent (Split-Path -Parent $MatlabExe)
    $toolboxDir = Join-Path $matlabRoot 'toolbox'

    # Mapping (display name) -> (cartella sotto toolbox\)
    $toolboxFolders = [ordered]@{
        'Mapping Toolbox' = 'map'
        'Image Processing Toolbox' = 'images'
        'Signal Processing Toolbox' = 'signal'
        'Statistics and Machine Learning Toolbox' = 'stats'
        'Parallel Computing Toolbox' = 'parallel'
    }

    $result = [ordered]@{}
    foreach ($displayName in $toolboxFolders.Keys) {
        $folder = $toolboxFolders[$displayName]
        $result[$displayName] = Test-Path -LiteralPath (Join-Path $toolboxDir $folder)
    }
    return $result
}

# Toolbox bloccanti per PHASE - se mancano, non si può fare PSI.
$Script:RequiredToolboxes = @('Mapping Toolbox', 'Image Processing Toolbox')

# SNAP: cerca gpt.exe in tutte le install standard.
# Returns absolute path to gpt.exe or $null.
function Find-Snap {
    $globs = @(
        'C:\Program Files\esa-snap*\bin\gpt.exe',
        'C:\Program Files\snap*\bin\gpt.exe',
        'C:\Program Files (x86)\esa-snap*\bin\gpt.exe',
        "$env:USERPROFILE\esa-snap*\bin\gpt.exe",
        "$env:LOCALAPPDATA\Programs\esa-snap*\bin\gpt.exe"
    )
    foreach ($g in $globs) {
        $hit = Get-ChildItem -Path $g -ErrorAction SilentlyContinue |
            Sort-Object FullName -Descending | Select-Object -First 1
        if ($hit) { return $hit.FullName }
    }
    $cmd = Get-Command gpt.exe -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }
    return $null
}

# Python 3.11+: prova `py -3.X` per X=20..11 poi `python` su PATH.
# Esclude esplicitamente lo stub Microsoft Store sotto \WindowsApps\.
# Returns @{ Exe = '...'; Version = '3.11.9' } or $null.
function Find-Python {
    $candidates = @()

    # Try py launcher with explicit version pins
    if (Get-Command py.exe -ErrorAction SilentlyContinue) {
        for ($minor = 20; $minor -ge $Script:PythonMinMinor; $minor--) {
            $candidates += "py -3.$minor"
        }
        $candidates += 'py -3'
    }

    # Fallback to plain python and python3 on PATH
    foreach ($name in @('python', 'python3')) {
        $cmd = Get-Command $name -ErrorAction SilentlyContinue
        if ($cmd -and $cmd.Source -notlike '*\WindowsApps\*') {
            $candidates += "`"$($cmd.Source)`""
        }
    }

    foreach ($candidate in $candidates) {
        try {
            $resolved = & cmd /c "$candidate -c `"import sys; print(sys.executable + '|' + '%d.%d.%d' % sys.version_info[:3])`" 2>nul"
            if (-not $resolved) { continue }
            $parts = $resolved.Trim().Split('|')
            if ($parts.Count -ne 2) { continue }
            $exe = $parts[0]
            $ver = $parts[1]
            if ($exe -like '*\WindowsApps\*') { continue }
            $verParts = $ver.Split('.')
            $major = [int]$verParts[0]
            $minor = [int]$verParts[1]
            if ($major -lt $Script:PythonMinMajor) { continue }
            if ($major -eq $Script:PythonMinMajor -and $minor -lt $Script:PythonMinMinor) { continue }
            return @{ Exe = $exe; Version = $ver }
        } catch {
            continue
        }
    }
    return $null
}

function Find-Git {
    $cmd = Get-Command git.exe -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }
    # Common install paths
    foreach ($p in @(
        'C:\Program Files\Git\bin\git.exe',
        'C:\Program Files (x86)\Git\bin\git.exe',
        "$env:LOCALAPPDATA\Programs\Git\bin\git.exe",
        # Portable Git installato dall'installer (se gia' eseguito una volta)
        "$env:LOCALAPPDATA\PHASE\portable-git\bin\git.exe"
    )) {
        if (Test-Path $p) { return $p }
    }
    return $null
}

# GMT (Generic Mapping Tools): cerca gmt.exe in Program Files glob + PATH
# + dir portable dell'installer.
# Returns absolute path a gmt.exe oppure $null.
function Find-Gmt {
    # Portable estratto dall'installer (preferito - controllato per primo)
    $portable = Join-Path $env:LOCALAPPDATA 'PHASE\gmt\bin\gmt.exe'
    if (Test-Path $portable) { return $portable }

    $globs = @(
        'C:\Program Files\GMT*\bin\gmt.exe',
        'C:\Program Files (x86)\GMT*\bin\gmt.exe',
        "$env:LOCALAPPDATA\Programs\GMT*\bin\gmt.exe"
    )
    foreach ($g in $globs) {
        $hit = Get-ChildItem -Path $g -ErrorAction SilentlyContinue |
            Sort-Object FullName -Descending | Select-Object -First 1
        if ($hit) { return $hit.FullName }
    }
    $cmd = Get-Command gmt.exe -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }
    return $null
}

# Scarica e ESTRAE la versione portable di GMT (~182 MB zip, ~700 MB
# estratto). Necessario per TRAIN con tropo_method='a_gacos'.
#
# Strategia ZIP portable (no installer NSIS):
#   - Niente prompt UAC (estrazione in %LOCALAPPDATA%\PHASE\gmt)
#   - Niente dialog "Failed to add to PATH" (l'installer NSIS lo mostra
#     anche in silent mode quando manca admin per HKLM)
#   - Pattern identico a Portable Git
#   - Idempotente (skip se gmt.exe gia' presente)
function Install-GmtSilent {
    param(
        [Parameter(Mandatory)] [scriptblock]$StatusCallback,
        [Parameter(Mandatory)] [scriptblock]$ProgressCallback
    )
    $existing = Find-Gmt
    if ($existing) {
        & $StatusCallback "GMT already present: $existing"
        # Anche se gia' presente, mi assicuro che <bin> sia su PATH user
        $existingBin = Split-Path -Parent $existing
        Add-DirToUserPath -Dir $existingBin -StatusCallback $StatusCallback
        return $existing
    }

    $destDir = Join-Path $env:LOCALAPPDATA 'PHASE\gmt'
    $gmtExe  = Join-Path $destDir 'bin\gmt.exe'

    & $StatusCallback 'Downloading GMT 6.6.0 portable (~182 MB, may take 2-5 min)...'
    $tmpZip = Join-Path $env:TEMP "gmt-portable_$(Get-Random).zip"
    try {
        Get-RemoteFile -Url $Script:GmtZipUrl -OutFile $tmpZip -ProgressCallback $ProgressCallback | Out-Null
        & $StatusCallback 'Extracting GMT (~700 MB on disk)...'

        if (-not (Test-Path $destDir)) {
            New-Item -ItemType Directory -Path $destDir -Force | Out-Null
        }
        # Lo zip GMT ha bin/, share/, lib/ al top-level (niente sub-cartella)
        Add-Type -AssemblyName System.IO.Compression.FileSystem -ErrorAction SilentlyContinue
        [System.IO.Compression.ZipFile]::ExtractToDirectory($tmpZip, $destDir)

        if (-not (Test-Path $gmtExe)) {
            throw "GMT extraction completed but gmt.exe not found at $gmtExe"
        }

        $gmtBin = Split-Path -Parent $gmtExe
        Add-DirToUserPath -Dir $gmtBin -StatusCallback $StatusCallback
        # Aggiorno anche il PATH della sessione corrente
        $env:Path = "$env:Path;$gmtBin"

        & $StatusCallback "[OK] GMT portable installed: $gmtExe"
        return $gmtExe
    } finally {
        Remove-Item $tmpZip -Force -ErrorAction SilentlyContinue
    }
}

# Helper: aggiunge una directory al PATH user-scope se non gia' presente.
# Usato per GMT bin e potenzialmente altre tool che devono essere su PATH.
function Add-DirToUserPath {
    param(
        [Parameter(Mandatory)] [string]$Dir,
        [Parameter(Mandatory)] [scriptblock]$StatusCallback
    )
    $userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
    if ($userPath -and ($userPath -split ';' | Where-Object { $_ -ieq $Dir })) {
        & $StatusCallback "$Dir already on user PATH"
        return
    }
    $newPath = if ($userPath) { "$userPath;$Dir" } else { $Dir }
    [Environment]::SetEnvironmentVariable('Path', $newPath, 'User')
    & $StatusCallback "[OK] $Dir aggiunto al PATH user"
}

# Scarica Portable Git da git-for-windows e lo estrae in una sotto-cartella
# %LOCALAPPDATA%\PHASE\portable-git. NON tocca l'installazione di sistema -
# resta una install isolata usata solo dall'installer e da PHASE.
#
# Il .7z.exe e' self-extracting: lanciato con `-o<dir> -y` estrae silent
# senza UAC ne' prompt (binari Inno Setup-style con flag standard 7-zip SFX).
#
# Non richiede credenziali git (i 3 repo sono pubblici, clone HTTPS anonimo).
function Install-PortableGit {
    param(
        [Parameter(Mandatory)] [scriptblock]$StatusCallback,
        [Parameter(Mandatory)] [scriptblock]$ProgressCallback
    )
    $destDir = Join-Path $env:LOCALAPPDATA 'PHASE\portable-git'
    $gitExe  = Join-Path $destDir 'bin\git.exe'

    if (Test-Path $gitExe) {
        & $StatusCallback "Portable Git already present at $destDir"
        return $gitExe
    }

    & $StatusCallback 'Download Portable Git (~50 MB)...'
    $tmpExe = Join-Path $env:TEMP "PortableGit-installer_$(Get-Random).exe"
    try {
        Get-RemoteFile -Url $Script:GitPortableUrl -OutFile $tmpExe -ProgressCallback $ProgressCallback | Out-Null
        & $StatusCallback 'Extracting Portable Git...'

        if (-not (Test-Path $destDir)) {
            New-Item -ItemType Directory -Path $destDir -Force | Out-Null
        }
        # 7-zip SFX flags: -o<dir> directory output, -y assume yes
        # NB: il comando va passato senza spazio tra -o e il path.
        $proc = Start-Process -FilePath $tmpExe `
            -ArgumentList "-o`"$destDir`"", '-y' `
            -Wait -PassThru -NoNewWindow
        if ($proc.ExitCode -ne 0 -and -not (Test-Path $gitExe)) {
            throw "Estrazione Portable Git fallita (exit $($proc.ExitCode))"
        }

        if (-not (Test-Path $gitExe)) {
            throw "Portable Git extracted but git.exe not found at $gitExe"
        }
        & $StatusCallback "[OK] Portable Git installed in $destDir"
        return $gitExe
    } finally {
        Remove-Item $tmpExe -Force -ErrorAction SilentlyContinue
    }
}

# -----------------------------------------------------------------------------
# Install / action helpers
# -----------------------------------------------------------------------------

# Download a file with progress to a target. Uses BITS when available (faster),
# falls back to Invoke-WebRequest. Returns the local path.
function Get-RemoteFile {
    param(
        [Parameter(Mandatory)] [string]$Url,
        [Parameter(Mandatory)] [string]$OutFile,
        [scriptblock]$ProgressCallback
    )
    $dir = Split-Path -Parent $OutFile
    if (-not (Test-Path $dir)) { New-Item -ItemType Directory -Path $dir -Force | Out-Null }

    $req = [System.Net.HttpWebRequest]::Create($Url)
    $req.UserAgent = 'PHASE-Installer/1.0'
    $resp = $req.GetResponse()
    $total = $resp.ContentLength
    $stream = $resp.GetResponseStream()
    $fileStream = [System.IO.File]::Create($OutFile)
    # Use a much larger buffer (1 MB) so the read loop yields the thread
    # less frequently and the I/O has time to flush. Smaller buffers
    # (81 KB) cause the loop to spin so fast that DoEvents pumping can't
    # keep up on large downloads (1 GB SNAP -> 12 800 iterations).
    $buffer = New-Object byte[] 1048576
    $read = 0
    $totalRead = 0
    $lastPct = -1
    try {
        while (($read = $stream.Read($buffer, 0, $buffer.Length)) -gt 0) {
            $fileStream.Write($buffer, 0, $read)
            $totalRead += $read
            if ($total -gt 0) {
                $pct = [Math]::Min(100, [int](($totalRead / $total) * 100))
                # Pump the WPF/WinForms message loop every iteration so the
                # main thread stays responsive (no "not responding" title
                # bar on multi-minute downloads). Cost: a few ms per MB,
                # negligible vs the network throughput.
                [System.Windows.Forms.Application]::DoEvents()
                if ($pct -ne $lastPct -and $ProgressCallback) {
                    & $ProgressCallback $pct
                    $lastPct = $pct
                }
            }
        }
    } finally {
        $fileStream.Close()
        $stream.Close()
        $resp.Close()
    }
    return $OutFile
}

# Silent install of CPython from python.org. Returns the new python.exe path.
function Install-PythonSilent {
    param(
        [Parameter(Mandatory)] [scriptblock]$ProgressCallback,
        [Parameter(Mandatory)] [scriptblock]$StatusCallback
    )
    & $StatusCallback 'Downloading...'
    $tmp = Join-Path $env:TEMP 'phase-python-installer.exe'
    Get-RemoteFile -Url $Script:PythonUrl -OutFile $tmp -ProgressCallback $ProgressCallback | Out-Null

    & $StatusCallback 'Installazione silenziosa Python 3.11.9...'
    # InstallAllUsers=1 -> sotto Program Files (richiede UAC)
    # InstallAllUsers=0 -> sotto AppData\Local\Programs\Python\Python311 (no UAC)
    # Scegliamo per-user per evitare il prompt UAC dentro al wizard.
    $args = @(
        '/quiet',
        'InstallAllUsers=0',
        'PrependPath=1',
        'Include_test=0',
        'Include_pip=1',
        'Include_launcher=1'
    )
    $proc = Start-Process -FilePath $tmp -ArgumentList $args -Wait -PassThru
    Remove-Item -Path $tmp -Force -ErrorAction SilentlyContinue
    if ($proc.ExitCode -ne 0) {
        throw "Python installer exited with code $($proc.ExitCode)"
    }

    # Aggiorna PATH della sessione corrente (l'installer lo setta nel registro
    # ma il processo corrente eredita il vecchio PATH)
    $userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
    $machinePath = [Environment]::GetEnvironmentVariable('Path', 'Machine')
    $env:Path = "$machinePath;$userPath"

    & $StatusCallback 'Post-install detection...'
    Start-Sleep -Seconds 2
    $found = Find-Python
    if (-not $found) {
        throw "Python installed but not detected on PATH after install. Restart the wizard."
    }
    return $found
}

function Install-PythonPackages {
    param(
        [Parameter(Mandatory)] [string]$PythonExe,
        [Parameter(Mandatory)] [scriptblock]$StatusCallback
    )
    # Pacchetti Python richiesti da PHASE:
    #   openpyxl    -> export/lettura Excel
    #   requests    -> download HTTP
    #   asf_search  -> ricerca e download dati SAR da ASF (Sentinel-1)
    #   shapely     -> geometrie / footprint AOI
    $packages = @('openpyxl', 'requests', 'asf_search', 'shapely')
    foreach ($pkg in $packages) {
        & $StatusCallback "pip install $pkg..."
        $proc = Start-Process -FilePath $PythonExe -ArgumentList '-m', 'pip', 'install', '--upgrade', $pkg -Wait -PassThru -NoNewWindow
        if ($proc.ExitCode -ne 0) {
            throw "pip install $pkg fallito con exit code $($proc.ExitCode)"
        }
    }
}

# Runs the SNAP installer (install4j) silently. If no local installer path
# is given, downloads it from ESA's CDN first.
#
# install4j silent flags:
#   -q          quiet, use default for every variable
#   -overwrite  overwrite existing files without prompt
#
# Single UAC prompt (Program Files write). If silent fails (exit != 0)
# falls back to the interactive install4j wizard.
function Invoke-SnapInstaller {
    param(
        [string]$InstallerPath,         # optional: local pre-staged installer
        [scriptblock]$StatusCallback,
        [scriptblock]$ProgressCallback,         # optional: download progress %
        [scriptblock]$InstallStartCallback      # optional: signalled when the
                                                # download finishes and the
                                                # install4j run begins. UI uses
                                                # this to switch the progress
                                                # bar to indeterminate marquee.
    )
    $statusCb = if ($StatusCallback) { $StatusCallback } else { { param($m) Write-Host $m } }

    # If no installer path was given (or the bundled file doesn't exist),
    # download SNAP from the ESA CDN. ~1 GB, expect 2-5 minutes on a fast
    # connection. The progress callback drives the SnapDownloadProgress
    # bar in the UI.
    $downloadedPath = $null
    if (-not $InstallerPath -or -not (Test-Path $InstallerPath)) {
        & $statusCb 'No bundled SNAP installer found - downloading from ESA (~1 GB, 2-5 min)...'
        $downloadedPath = Join-Path $env:TEMP "esa-snap-installer_$(Get-Random).exe"
        try {
            Get-RemoteFile -Url $Script:SnapInstallerUrl -OutFile $downloadedPath -ProgressCallback $ProgressCallback | Out-Null
            $InstallerPath = $downloadedPath
            & $statusCb "Downloaded $([math]::Round((Get-Item $InstallerPath).Length / 1MB)) MB from ESA"
        } catch {
            if ($downloadedPath -and (Test-Path $downloadedPath)) {
                Remove-Item $downloadedPath -Force -ErrorAction SilentlyContinue
            }
            throw "Unable to download SNAP installer from $($Script:SnapInstallerUrl): $($_.Exception.Message)"
        }
    }

    try {
        & $statusCb 'Running SNAP installer silently (~3-5 min, UAC required)...'
        if ($InstallStartCallback) { & $InstallStartCallback }

        # Attempt 1: silent install
        $proc = Start-Process -FilePath $InstallerPath `
            -ArgumentList '-q', '-overwrite' `
            -Wait -PassThru
        if ($proc.ExitCode -eq 0) {
            & $statusCb 'SNAP silent install completed (exit 0)'
            return 0
        }

        & $statusCb "SNAP silent install returned exit code $($proc.ExitCode) - re-launching in interactive mode"
        # Attempt 2: fallback to interactive wizard
        $proc = Start-Process -FilePath $InstallerPath -Wait -PassThru
        return $proc.ExitCode
    } finally {
        # Clean up the downloaded installer (only if we created it)
        if ($downloadedPath -and (Test-Path $downloadedPath)) {
            Remove-Item $downloadedPath -Force -ErrorAction SilentlyContinue
        }
    }
}

# git clone con progress callback (parse output di --progress).
function Invoke-GitClone {
    param(
        [Parameter(Mandatory)] [string]$GitExe,
        [Parameter(Mandatory)] [string]$Repo,
        [Parameter(Mandatory)] [string]$Branch,
        [Parameter(Mandatory)] [string]$Destination,
        [Parameter(Mandatory)] [scriptblock]$StatusCallback
    )
    if (Test-Path (Join-Path $Destination '.git')) {
        # Repo gia' presente: forza l'allineamento al branch corretto del
        # remote configurato. Usato sia per riprese di install interrotte
        # sia per branch rinominati upstream (es. windows-port/main -> main).
        # Strategia:
        #   1. remote set-url origin <Repo>          (gestisce fork swap)
        #   2. reset --hard HEAD                     (scarta modifiche locali
        #                                             tracked - tipicamente
        #                                             le patch .mlapp precedenti)
        #   3. fetch origin <Branch>
        #   4. checkout -B <Branch> FETCH_HEAD       (anche su storia divergente
        #                                             - bypassa il ff-only)
        & $StatusCallback "Repo already present at $Destination - updating to origin/$Branch..."
        $null = Start-Process -FilePath $GitExe -ArgumentList @('-C', $Destination, 'remote', 'set-url', 'origin', $Repo) -Wait -PassThru -NoNewWindow
        $null = Start-Process -FilePath $GitExe -ArgumentList @('-C', $Destination, 'reset', '--hard', 'HEAD') -Wait -PassThru -NoNewWindow
        $pFetch = Start-Process -FilePath $GitExe -ArgumentList @('-C', $Destination, 'fetch', 'origin', $Branch) -Wait -PassThru -NoNewWindow
        if ($pFetch.ExitCode -ne 0) {
            & $StatusCallback "git fetch fallito (exit $($pFetch.ExitCode)) - mantengo lo stato attuale"
            return
        }
        $pCheckout = Start-Process -FilePath $GitExe -ArgumentList @('-C', $Destination, 'checkout', '-B', $Branch, 'FETCH_HEAD') -Wait -PassThru -NoNewWindow
        if ($pCheckout.ExitCode -eq 0) {
            & $StatusCallback "Repo allineato al branch $Branch (HEAD da remote)"
        } else {
            & $StatusCallback "git checkout fallito (exit $($pCheckout.ExitCode))"
        }
        return
    }
    if (Test-Path $Destination) {
        throw "Cartella $Destination esiste ma non e' un repo git. Spostala o cancellala manualmente."
    }

    & $StatusCallback "git clone $Repo (branch $Branch)..."
    $proc = Start-Process -FilePath $GitExe `
        -ArgumentList 'clone', '--branch', $Branch, '--single-branch', $Repo, $Destination `
        -Wait -PassThru -NoNewWindow
    if ($proc.ExitCode -ne 0) {
        throw "git clone $Repo fallito con exit code $($proc.ExitCode)"
    }
}

# Scrive %APPDATA%\PHASE\python.txt con il path al Python interpreter.
# Letto da StaMPS\bin\mt_prep_snap.bat:27 per bypassare la ricerca di `py -3`.
function Set-PhasePythonConfig {
    param([Parameter(Mandatory)] [string]$PythonExe)
    $dir = Join-Path $env:APPDATA 'PHASE'
    if (-not (Test-Path $dir)) { New-Item -ItemType Directory -Path $dir -Force | Out-Null }
    $file = Join-Path $dir 'python.txt'
    Set-Content -Path $file -Value $PythonExe -Encoding ASCII -NoNewline
}

# setx MATLAB_EXE (user scope, persistente).
function Set-MatlabEnvVar {
    param([Parameter(Mandatory)] [string]$MatlabExe)
    [Environment]::SetEnvironmentVariable('MATLAB_EXE', $MatlabExe, 'User')
}

# Lancia MATLAB in batch per:
#   1. addpath + savepath su StaMPS\matlab, StaMPS\matlab_compat e TRAIN\matlab
#   2. (PhaseRoot) generare input_StaMPS.mat precompilato con installation_folder
#      e project_path gia' settati, cosi' l'utente apre il .mlapp e clicca
#      Load senza dover sfogliare i due path manualmente.
# Fallback graceful: se MATLAB rifiuta -batch (vecchie versioni) o se la
# licenza non è ancora attivata, scrive un warning nel log ma non blocca.
function Invoke-MatlabSavePath {
    param(
        [Parameter(Mandatory)] [string]$MatlabExe,
        [Parameter(Mandatory)] [string]$StampsRoot,
        [string]$TrainRoot,
        [string]$PhaseRoot,
        [string]$PythonExe,
        [string]$SnapGpt,
        [Parameter(Mandatory)] [scriptblock]$StatusCallback
    )
    # Costruisco lo script MATLAB come file .m temporaneo. Più affidabile
    # che passare statement inline a -batch (evita escaping di virgolette
    # e apostrofi attraverso PowerShell -> cmd -> matlab.exe).
    $stampsM = $StampsRoot.Replace('\','/')
    $mLines = @(
        "try"
        "    addpath(genpath('$stampsM/matlab'));"
        "    addpath(genpath('$stampsM/matlab_compat'));"
    )
    if ($TrainRoot) {
        $trainM = $TrainRoot.Replace('\','/')
        $mLines += "    addpath(genpath('$trainM/matlab'));"
    }
    $mLines += @(
        "    savepath;"
    )

    # Step 2: scrittura input_StaMPS.mat precompilato (60+ var con i default
    # del .mlapp + i due path dell'installer). Le 60+ variabili sono
    # necessarie perche' StartButtonPushed (document.xml:623) fa load di
    # *tutte* le variabili - se ne manca una, il load esplode.
    if ($PhaseRoot) {
        $phaseM = $PhaseRoot.Replace('\','/')
        $matFile = "$phaseM/PHASE_Preprocessing/input_StaMPS.mat"
        $installFolder = "$phaseM/StaMPS"
        # NB: project_path in PHASE_StaMPS e' la PHASE root, non la sottocartella
        # PHASE_Preprocessing. PHASE_StaMPS appende lui stesso "/PHASE_Preprocessing/INSAR_<date>"
        # (vedi document.xml:682). Scrivere "$phaseM/PHASE_Preprocessing" qui
        # duplicava il segmento => path inesistente.
        $projectPath = "$phaseM"
        $mLines += @(
            "    stamps_preparation = 0;"
            "    installation_folder = '$installFolder';"
            "    project_path = '$projectPath';"
            "    amplitude_threshold = 0.40;"
            "    master_date = '20200305';"
            "    export_name = 'filename';"
            "    time_span = 0;"
            "    year_0 = 2019;"
            "    month_0 = 7;"
            "    day_0 = 9;"
            "    utc_time = '10:30';"
            "    train_flag = 0;"
            "    stamps_first_step = '1';"
            "    stamps_last_step = '7';"
            "    n_cores = 4;"
            "    heading = 346.18;"
            "    lambda = 0.055465763;"
            "    max_topo_err = 16;"
            "    filter_grid_size = 40;"
            "    filter_weighting = 'P-square';"
            "    gamma_max_iterations = 7;"
            "    gamma_change_convergence = 0.0050;"
            "    gamma_stdev_reject = 0;"
            "    quick_est_gamma_flag = 'y';"
            "    small_baseline_flag = 'n';"
            "    clap_win = 16;"
            "    clap_alpha = 1;"
            "    clap_beta = 0.3000;"
            "    clap_low_pass_wavelength = 800;"
            "    select_method = 'PERCENT';"
            "    percent_rand = 1;"
            "    weed_standard_dev = 1;"
            "    weed_neighbours = 'y';"
            "    weed_zero_elevation = 'n';"
            "    weed_max_noise = Inf;"
            "    merge_resample_size = 0;"
            "    merge_standard_dev = Inf;"
            "    unwrap_grid_size = 20;"
            "    unwrap_gold_n_win = 16;"
            "    unwrap_method = '3D';"
            "    unwrap_gold_alpha = 0.8;"
            "    unwrap_alpha = 8;"
            "    unwrap_spatial_cost_func_flag = 'n';"
            "    unwrap_prefilter_flag = 'y';"
            "    unwrap_patch_phase = 'n';"
            "    unwrap_la_error_flag = 'y';"
            "    unwrap_hold_good_values = 'y';"
            # Default a 'n': l'utente puo' attivare la correzione troposferica
            # in GUI dopo aver installato TRAIN + i dati meteo (NCEP / ERA5 / GACOS).
            # Con 'y' di default, step 7 (ps_calc_scla) finisce in keyboard
            # mode su un'installazione fresca senza TRAIN ancora configurato.
            "    subtr_tropo = 'n';"
            "    tropo_method = 'a_linear';"
            "    select_reest_gamma_flag = 'y';"
            "    drop_ifg_index = '[]';"
            "    scla_deramp = 'y';"
            "    scla_method = 'L2';"
            "    scla_drop_index = '[]';"
            "    scn_wavelength = 50;"
            "    scn_kriging_flag = 'n';"
            "    ref_centre_lonlat = [0.0 0.0];"
            "    ref_radius = 0;"
            "    ref_velocity = 0;"
            "    plot_s = 15;"
            "    ref_centre_lonlat_w = [0.0 0.0];"
            "    ref_radius_w = 0;"
            "    ph_output = 'unwrapped';"
            "    save('$matFile', 'stamps_preparation', 'installation_folder', 'project_path', 'amplitude_threshold', 'master_date', 'export_name', 'time_span', 'year_0', 'month_0', 'day_0', 'utc_time', 'train_flag', 'stamps_first_step', 'stamps_last_step', 'n_cores', 'heading', 'lambda', 'max_topo_err', 'filter_grid_size', 'filter_weighting', 'gamma_max_iterations', 'gamma_change_convergence', 'gamma_stdev_reject', 'quick_est_gamma_flag', 'small_baseline_flag', 'clap_win', 'clap_alpha', 'clap_beta', 'clap_low_pass_wavelength', 'select_method', 'percent_rand', 'weed_standard_dev', 'weed_neighbours', 'weed_zero_elevation', 'weed_max_noise', 'merge_resample_size', 'merge_standard_dev', 'unwrap_grid_size', 'unwrap_gold_n_win', 'unwrap_method', 'unwrap_gold_alpha', 'unwrap_alpha', 'unwrap_spatial_cost_func_flag', 'unwrap_prefilter_flag', 'unwrap_patch_phase', 'unwrap_la_error_flag', 'unwrap_hold_good_values', 'subtr_tropo', 'tropo_method', 'select_reest_gamma_flag', 'drop_ifg_index', 'scla_deramp', 'scla_method', 'scla_drop_index', 'scn_wavelength', 'scn_kriging_flag', 'ref_centre_lonlat', 'ref_radius', 'ref_velocity', 'plot_s', 'ref_centre_lonlat_w', 'ref_radius_w', 'ph_output', '-mat');"
        )
    }

    # Step 3: scrittura input_preprocessing.mat precompilato (24 var; configurazione
    # SEN di default - costellazione piu' usata. Le 2 var di installazione sono
    # python (=PythonExe dell'installer) e gptbin_path (=SnapGpt). Le altre 22
    # sono i default UI del .mlapp - l'utente le modifica per ogni dataset.
    # NB: il LoadButton di PHASE_Preprocessing.mlapp legge solo le var presenti,
    # quindi se l'utente cambia a CSK ricarica le 20 var CSK dal save successivo.
    if ($PhaseRoot -and $PythonExe -and $SnapGpt) {
        $phaseM = $PhaseRoot.Replace('\','/')
        $prepMatFile = "$phaseM/PHASE_Preprocessing/input_preprocessing.mat"
        $pyM = $PythonExe.Replace('\','/')
        $gptM = $SnapGpt.Replace('\','/')
        $mLines += @(
            "    constellation = 'SEN';"
            "    python = '$pyM';"
            "    gptbin_path = '$gptM';"
            "    images_download = 1;"
            "    master_date = '20200722';"
            "    auto_master = 1;"
            "    master_processing = 0;"
            "    polarisation = 'VV';"
            "    lon_min = -180.000;"
            "    lat_min = -90.000;"
            "    lon_max = 180.000;"
            "    lat_max = 90.000;"
            "    slaves_removal = 1;"
            "    dem_name = 'SRTM 1Sec HGT';"
            "    dem_file = '';"
            "    dem_name_coreg = 'SRTM 1Sec HGT';"
            "    dem_file_coreg = '';"
            "    dem_resampling = 'NEAREST_NEIGHBOUR';"
            "    first_step = 1;"
            "    coherence_tc = 0;"
            "    epsg_code = 32633;"
            "    cpu = 8;"
            "    cache = '26G';"
            "    save('$prepMatFile', 'constellation', 'python', 'images_download', 'master_date', 'auto_master', 'master_processing', 'polarisation', 'lon_min', 'lat_min', 'lon_max', 'lat_max', 'slaves_removal', 'dem_name', 'dem_file', 'dem_name_coreg', 'dem_file_coreg', 'dem_resampling', 'first_step', 'coherence_tc', 'epsg_code', 'gptbin_path', 'cpu', 'cache', '-mat');"
        )
    }

    $mLines += @(
        "    fprintf('PHASE_INSTALLER_SAVEPATH_OK\n');"
        "catch err"
        "    fprintf('PHASE_INSTALLER_SAVEPATH_ERR: %s\n', err.message);"
        "end"
        "exit;"
    )
    $tmpScript = Join-Path $env:TEMP "phase_savepath_$(Get-Random).m"
    Set-Content -Path $tmpScript -Value ($mLines -join "`r`n") -Encoding ASCII

    & $StatusCallback 'MATLAB savepath + writing input_StaMPS.mat (may take 30-60s)...'

    # Uso System.Diagnostics.Process direttamente. Start-Process con
    # -NoNewWindow + -Wait + -RedirectStandardOutput non funziona affidabile
    # con app GUI come matlab.exe (il process parent puo' exit prima del
    # batch completion).
    $psi = New-Object System.Diagnostics.ProcessStartInfo
    $psi.FileName = $MatlabExe
    $psi.Arguments = "-batch ""run('$($tmpScript.Replace('\','/'))')"""
    $psi.UseShellExecute = $false
    $psi.RedirectStandardOutput = $true
    $psi.RedirectStandardError = $true
    $psi.CreateNoWindow = $true

    $proc = New-Object System.Diagnostics.Process
    $proc.StartInfo = $psi
    [void]$proc.Start()
    $stdout = $proc.StandardOutput.ReadToEnd()
    $stderr = $proc.StandardError.ReadToEnd()
    $proc.WaitForExit()
    $exitCode = $proc.ExitCode

    Remove-Item $tmpScript -ErrorAction SilentlyContinue

    $combined = $stdout + "`n" + $stderr
    $success = ($combined -match 'PHASE_INSTALLER_SAVEPATH_OK')
    return @{
        Success = $success
        Output  = $combined.Trim()
        ExitCode = $exitCode
    }
}

# Modifica un .mlapp clonato per aggiungere auto-load nello startupFcn.
# Quando l'utente apre il .mlapp con doppio click, se input_*.mat esiste
# l'app chiama in automatico il callback del LoadButton e tutti i campi
# si popolano (compresi installation_folder e project_path) senza che
# l'utente debba andare al tab Save/Load e cliccare Load.
#
# I .mlapp sono archivi zip con struttura OPC: i path interni usano '/'
# come separator (non '\') e l'ordine delle entry e' significativo. Quindi
# qui re-zippiamo manualmente preservando ordine + separator originale.
function Invoke-MlappAutoLoadPatch {
    param(
        [Parameter(Mandatory)] [string]$MlappPath,
        [string]$MatFileRelative,
        [string]$Anchor = 'cd(currentFolder);',
        [string]$InjectBlock,
        [Parameter(Mandatory)] [scriptblock]$StatusCallback
    )
    if (-not (Test-Path -LiteralPath $MlappPath)) {
        & $StatusCallback "Mlapp not found: $MlappPath (skipping)"
        return $false
    }

    Add-Type -AssemblyName System.IO.Compression -ErrorAction SilentlyContinue
    Add-Type -AssemblyName System.IO.Compression.FileSystem -ErrorAction SilentlyContinue

    $extractTmp = Join-Path $env:TEMP "phase_mlapp_patch_$(Get-Random)"
    New-Item -ItemType Directory -Path $extractTmp -Force | Out-Null

    try {
        # Salva ordine originale delle entry (cruciale per OPC/MATLAB)
        $orderedEntries = @()
        $zipRead = [System.IO.Compression.ZipFile]::OpenRead($MlappPath)
        $orderedEntries = $zipRead.Entries | ForEach-Object { $_.FullName }
        $zipRead.Dispose()

        # Estrai tutto in tmp
        [System.IO.Compression.ZipFile]::ExtractToDirectory($MlappPath, $extractTmp)

        # Modifica matlab/document.xml
        $docXml = Join-Path $extractTmp 'matlab\document.xml'
        $content = [System.IO.File]::ReadAllText($docXml, [System.Text.Encoding]::UTF8)

        # Verifica che la patch non sia gia' stata applicata (idempotenza)
        if ($content.Contains('AUTO-LOAD (PHASE installer)')) {
            & $StatusCallback "$([System.IO.Path]::GetFileName($MlappPath)): patch already applied, skipping"
            return $true
        }

        # Anchor: prima occorrenza della stringa anchor
                $idx = $content.IndexOf($Anchor)
        if ($idx -lt 0) {
            & $StatusCallback "Warning: Anchor '$Anchor' not found. Skipping patch auto-load."
            return $false
        }
        $insertPoint = $idx + $Anchor.Length

        # InjectBlock: o quello passato dal chiamante, o il default auto-load
        # standard che chiama LoadButton se il .mat esiste.
        if ($InjectBlock) {
            $blockToInsert = $InjectBlock
        } else {
            $blockToInsert = @"


            % AUTO-LOAD (PHASE installer): applica i default precompilati se il
            % .mat di default esiste. L'utente non deve cliccare manualmente Load.
            try
                if exist('$MatFileRelative', 'file') == 2 && isprop(app, 'LoadButton')
                    feval(app.LoadButton.ButtonPushedFcn, app.LoadButton, struct());
                end
            catch
                % non blocca lo startup se il load fallisce
            end
"@
        }

        $patched = $content.Substring(0, $insertPoint) + $blockToInsert + $content.Substring($insertPoint)
        # CRITICAL: scrivere UTF-8 SENZA BOM. App Designer di MATLAB rifiuta
        # i .mlapp con BOM nel document.xml (la classe viene parsata ma TUTTE
        # le callback risultano "non definite" - bug ostico da diagnosticare).
        # System.Text.Encoding.UTF8 default include BOM; usiamo new($false).
        $utf8NoBom = New-Object System.Text.UTF8Encoding $false
        [System.IO.File]::WriteAllText($docXml, $patched, $utf8NoBom)

        # Re-zippa preservando ordine + separator '/'
        Remove-Item -LiteralPath $MlappPath -Force
        $zip = [System.IO.Compression.ZipFile]::Open($MlappPath, 'Create')
        try {
            foreach ($name in $orderedEntries) {
                $fsPath = Join-Path $extractTmp ($name -replace '/', '\')
                if (-not (Test-Path -LiteralPath $fsPath)) {
                    & $StatusCallback "Warning: entry '$name' missing during repack"
                    continue
                }
                $entry = $zip.CreateEntry($name, [System.IO.Compression.CompressionLevel]::Optimal)
                $stream = $entry.Open()
                $bytes = [System.IO.File]::ReadAllBytes($fsPath)
                $stream.Write($bytes, 0, $bytes.Length)
                $stream.Close()
            }
        } finally {
            $zip.Dispose()
        }

        & $StatusCallback "$([System.IO.Path]::GetFileName($MlappPath)): auto-load aggiunto (anchor='$MatFileRelative')"
        return $true

    } finally {
        Remove-Item -Recurse -Force $extractTmp -ErrorAction SilentlyContinue
    }
}

# Genera un project.conf template con GPTBIN_PATH e GRAPHSFOLDER preconfigurati.
function Write-ProjectConfTemplate {
    param(
        [Parameter(Mandatory)] [string]$InstallDir,
        [Parameter(Mandatory)] [string]$SnapGpt
    )
    $template = @"
######### CONFIGURATION FILE per snap2stamps + PHASE - generato dall'installer ##########
# Aggiorna PROJECTFOLDER e MASTER per ogni nuovo dataset.

PROJECTFOLDER = $($InstallDir.Replace('\','/'))/PHASE_Preprocessing
GRAPHSFOLDER  = $($InstallDir.Replace('\','/'))/PHASE_Preprocessing/snap2stamps/graphs
GPTBIN_PATH   = $($SnapGpt.Replace('\','/'))

CACHE = 8G
CPU = 4

# Master: path al .dim splittato del master (compilato a runtime dall'app)
MASTER =

# AOI bounding box (richiesto dal pre-cache SRTM 3Sec)
LONMIN = 0.0
LATMIN = 0.0
LONMAX = 0.0
LATMAX = 0.0

SWATHS = IW1,IW2,IW3
POLARISATION = VV

DEMNAME = Copernicus 30m Global DEM
DEMFILE =
DEMRESAMPLING = BICUBIC_INTERPOLATION
"@
    $confPath = Join-Path $InstallDir 'project.conf.template'
    Set-Content -Path $confPath -Value $template -Encoding UTF8
}

# Scarica i 9 .exe di StaMPS dalla release "windows-port-bins-v1" di
# pyccino/StaMPS (asset stamps-win64-binaries.zip, ~4 MB compressed):
#   - 7 core tools (calamp, cpxsum, pscphase, pscdem, psclonlat, selpsc_patch,
#     selsbc_patch) → StaMPS/bin/
#   - triangle.exe → StaMPS/external/triangle/bin/
#   - snaphu.exe   → StaMPS/external/snaphu/bin/
# Necessari per il workflow PSI (mt_prep_snap, ps_load_initial_gamma) e per
# l'unwrap statistical-cost (ps_unwrap → uw_3d → uw_stat_costs invoca snaphu).
# Senza, StaMPS non puo' processare nulla (mt_prep fallisce al primo step,
# l'unwrap senza snaphu cadrebbe sul vecchio uw_nosnaphu che entra in loop
# infinito su AOI sparse).
#
# Questo e' il comportamento unico e di default: non tentiamo piu' di
# compilare Triangle/snaphu da sorgente con install-windows.ps1 (richiedeva
# MSVC/CMake e in pratica falliva sempre, attivando comunque questo fallback).
function Invoke-StampsBinariesDownload {
    param(
        [Parameter(Mandatory)] [string]$StampsRoot,
        [Parameter(Mandatory)] [scriptblock]$StatusCallback
    )
    $binDir = Join-Path $StampsRoot 'bin'
    $triangleBinDir = Join-Path $StampsRoot 'external\triangle\bin'
    $snaphuBinDir = Join-Path $StampsRoot 'external\snaphu\bin'
    if (-not (Test-Path $binDir)) {
        New-Item -ItemType Directory -Path $binDir -Force | Out-Null
    }
    if (-not (Test-Path $triangleBinDir)) {
        New-Item -ItemType Directory -Path $triangleBinDir -Force | Out-Null
    }
    if (-not (Test-Path $snaphuBinDir)) {
        New-Item -ItemType Directory -Path $snaphuBinDir -Force | Out-Null
    }

    # Mappa dei .exe richiesti al loro dest path. Tutti i .exe vivono al
    # root della zip; il routing per-nome decide se vanno in StaMPS/bin/
    # (tools StaMPS) o in StaMPS/external/triangle/bin/ (Triangle e' una
    # dipendenza esterna referenziata cosi' da mt_prep_snap.bat e
    # StaMPS_CONFIG.ps1).
    $exeDest = @{
        'calamp.exe'       = $binDir
        'cpxsum.exe'       = $binDir
        'pscphase.exe'     = $binDir
        'pscdem.exe'       = $binDir
        'psclonlat.exe'    = $binDir
        'selpsc_patch.exe' = $binDir
        'selsbc_patch.exe' = $binDir
        'triangle.exe'     = $triangleBinDir
        'snaphu.exe'       = $snaphuBinDir
    }
    $required = $exeDest.Keys
    $missing = $required | Where-Object { -not (Test-Path (Join-Path $exeDest[$_] $_)) }
    if ($missing.Count -eq 0) {
        & $StatusCallback "All $($required.Count) StaMPS .exe already present, skipping download"
        return $true
    }
    & $StatusCallback "Mancano $($missing.Count)/$($required.Count) .exe StaMPS, scarico stamps-win64-binaries.zip..."

    $url = 'https://github.com/pyccino/StaMPS/releases/download/windows-port-bins-v1/stamps-win64-binaries.zip'
    $tmpZip = Join-Path $env:TEMP "stamps-win64-binaries_$(Get-Random).zip"

    try {
        # Download
        $req = [System.Net.HttpWebRequest]::Create($url)
        $req.UserAgent = 'PHASE-Installer/1.1'
        $resp = $req.GetResponse()
        $total = $resp.ContentLength
        $in = $resp.GetResponseStream()
        $out = [System.IO.File]::Create($tmpZip)
        $buf = New-Object byte[] 81920
        $totalRead = 0
        while (($read = $in.Read($buf, 0, $buf.Length)) -gt 0) {
            $out.Write($buf, 0, $read)
            $totalRead += $read
        }
        $out.Close(); $in.Close(); $resp.Close()
        & $StatusCallback "Download complete ($([math]::Round($totalRead/1MB,1)) MB)"

        # Estrai i .exe nella loro dest dir (StaMPS/bin/ o
        # StaMPS/external/triangle/bin/ a seconda del nome).
        Add-Type -AssemblyName System.IO.Compression.FileSystem -ErrorAction SilentlyContinue
        $zip = [System.IO.Compression.ZipFile]::OpenRead($tmpZip)
        try {
            foreach ($entry in $zip.Entries) {
                if ($exeDest.ContainsKey($entry.Name)) {
                    $dest = Join-Path $exeDest[$entry.Name] $entry.Name
                    # Sovrascrive (i .exe possono evolvere tra release)
                    if (Test-Path $dest) { Remove-Item -LiteralPath $dest -Force }
                    [System.IO.Compression.ZipFileExtensions]::ExtractToFile($entry, $dest)
                }
            }
        } finally {
            $zip.Dispose()
        }

        # Verifica finale
        $stillMissing = $required | Where-Object { -not (Test-Path (Join-Path $exeDest[$_] $_)) }
        if ($stillMissing.Count -eq 0) {
            & $StatusCallback "$($required.Count)/$($required.Count) .exe StaMPS estratti (bin/ + external/triangle/bin/ + external/snaphu/bin/)"
            return $true
        } else {
            & $StatusCallback "ERRORE: dopo l'estrazione mancano ancora: $($stillMissing -join ', ')"
            return $false
        }
    } catch {
        & $StatusCallback "ERRORE download/estrazione binari StaMPS: $($_.Exception.Message)"
        return $false
    } finally {
        Remove-Item $tmpZip -Force -ErrorAction SilentlyContinue
    }
}

# -----------------------------------------------------------------------------
# WPF Wizard
# -----------------------------------------------------------------------------

[xml]$xaml = @'
<Window xmlns="http://schemas.microsoft.com/winfx/2006/xaml/presentation"
        xmlns:x="http://schemas.microsoft.com/winfx/2006/xaml"
        Title="PHASE Installer"
        Width="980" Height="700"
        WindowStartupLocation="CenterScreen"
        ResizeMode="NoResize"
        AllowsTransparency="False"
        Background="#FFFCFCFE"
        FontFamily="Segoe UI Variable, Segoe UI">
    <Window.Resources>
        <!-- ============== PALETTE PHASE — light scientific ==============
             surface:  #FCFCFE  base / window (off-white)
             elevated: #FFFFFF  cards, textboxes (pure white)
             sidebar:  #F4F6FB  sidebar tinted off-white (separazione)
             border:   #E4E8F0  hairline default
             borderHi: #C8D0E0  hairline hover
             text:     #0F1430  primary text (deep ink navy)
             text2:    #4A5168  secondary text
             text3:    #8C95B8  tertiary / placeholder / muted
             accent:   #1A4FE0  primary blue (vivido, leggibile su bianco)
             accent2:  #D1397E  magenta del logo
             accentLt: #E8EEFF  hover/selected very light blue
             success:  #2DBA6E
             warning:  #E89C2D  -->

        <!-- Solid blue primary button con glow soft (Avanti / Avvia) -->
        <Style x:Key="PrimaryButton" TargetType="Button">
            <Setter Property="Foreground" Value="White"/>
            <Setter Property="Background" Value="#1A4FE0"/>
            <Setter Property="BorderThickness" Value="0"/>
            <Setter Property="Padding" Value="26,11"/>
            <Setter Property="MinWidth" Value="124"/>
            <Setter Property="FontSize" Value="12"/>
            <Setter Property="FontWeight" Value="Bold"/>
            <Setter Property="Margin" Value="6,0"/>
            <Setter Property="Cursor" Value="Hand"/>
            <Setter Property="Template">
                <Setter.Value>
                    <ControlTemplate TargetType="Button">
                        <Border x:Name="border" Background="{TemplateBinding Background}" CornerRadius="4" Padding="{TemplateBinding Padding}">
                            <Border.Effect>
                                <DropShadowEffect Color="#FF1A4FE0" BlurRadius="18" ShadowDepth="0" Opacity="0.35"/>
                            </Border.Effect>
                            <ContentPresenter HorizontalAlignment="Center" VerticalAlignment="Center"/>
                        </Border>
                        <ControlTemplate.Triggers>
                            <Trigger Property="IsMouseOver" Value="True">
                                <Setter TargetName="border" Property="Background" Value="#2B62F0"/>
                                <Setter TargetName="border" Property="Effect">
                                    <Setter.Value>
                                        <DropShadowEffect Color="#FF1A4FE0" BlurRadius="30" ShadowDepth="0" Opacity="0.55"/>
                                    </Setter.Value>
                                </Setter>
                            </Trigger>
                            <Trigger Property="IsPressed" Value="True">
                                <Setter TargetName="border" Property="Background" Value="#143ECC"/>
                            </Trigger>
                            <Trigger Property="IsEnabled" Value="False">
                                <Setter TargetName="border" Property="Background" Value="#D6DCE8"/>
                                <Setter Property="Foreground" Value="#8C95B8"/>
                                <Setter TargetName="border" Property="Effect" Value="{x:Null}"/>
                            </Trigger>
                        </ControlTemplate.Triggers>
                    </ControlTemplate>
                </Setter.Value>
            </Setter>
        </Style>

        <!-- Ghost outlined secondary button (Indietro / Sfoglia / Annulla) -->
        <Style x:Key="SecondaryButton" TargetType="Button">
            <Setter Property="Background" Value="White"/>
            <Setter Property="Foreground" Value="#0F1430"/>
            <Setter Property="BorderBrush" Value="#D6DCE8"/>
            <Setter Property="BorderThickness" Value="1"/>
            <Setter Property="Padding" Value="22,10"/>
            <Setter Property="MinWidth" Value="110"/>
            <Setter Property="FontSize" Value="12"/>
            <Setter Property="FontWeight" Value="SemiBold"/>
            <Setter Property="Margin" Value="6,0"/>
            <Setter Property="Cursor" Value="Hand"/>
            <Setter Property="Template">
                <Setter.Value>
                    <ControlTemplate TargetType="Button">
                        <Border x:Name="border" Background="{TemplateBinding Background}" BorderBrush="{TemplateBinding BorderBrush}" BorderThickness="{TemplateBinding BorderThickness}" CornerRadius="4" Padding="{TemplateBinding Padding}">
                            <ContentPresenter HorizontalAlignment="Center" VerticalAlignment="Center"/>
                        </Border>
                        <ControlTemplate.Triggers>
                            <Trigger Property="IsMouseOver" Value="True">
                                <Setter TargetName="border" Property="Background" Value="#F4F6FC"/>
                                <Setter TargetName="border" Property="BorderBrush" Value="#1A4FE0"/>
                                <Setter Property="Foreground" Value="#1A4FE0"/>
                            </Trigger>
                            <Trigger Property="IsEnabled" Value="False">
                                <Setter Property="Foreground" Value="#B0B8C8"/>
                                <Setter TargetName="border" Property="BorderBrush" Value="#EEF1F8"/>
                            </Trigger>
                        </ControlTemplate.Triggers>
                    </ControlTemplate>
                </Setter.Value>
            </Setter>
        </Style>

        <Style TargetType="Button" BasedOn="{StaticResource SecondaryButton}"/>

        <Style TargetType="TextBlock">
            <Setter Property="FontFamily" Value="Segoe UI Variable, Segoe UI"/>
            <Setter Property="Foreground" Value="#0F1430"/>
        </Style>

        <Style TargetType="TextBox">
            <Setter Property="Padding" Value="14,11"/>
            <Setter Property="FontFamily" Value="JetBrains Mono, Cascadia Code, Consolas"/>
            <Setter Property="FontSize" Value="12"/>
            <Setter Property="Foreground" Value="#0F1430"/>
            <Setter Property="CaretBrush" Value="#1A4FE0"/>
            <Setter Property="SelectionBrush" Value="#1A4FE0"/>
            <Setter Property="BorderBrush" Value="#E4E8F0"/>
            <Setter Property="BorderThickness" Value="1"/>
            <Setter Property="Template">
                <Setter.Value>
                    <ControlTemplate TargetType="TextBox">
                        <Border x:Name="border" Background="White" BorderBrush="{TemplateBinding BorderBrush}" BorderThickness="{TemplateBinding BorderThickness}" CornerRadius="3">
                            <ScrollViewer x:Name="PART_ContentHost" Padding="{TemplateBinding Padding}"/>
                        </Border>
                        <ControlTemplate.Triggers>
                            <Trigger Property="IsFocused" Value="True">
                                <Setter TargetName="border" Property="BorderBrush" Value="#1A4FE0"/>
                            </Trigger>
                            <Trigger Property="IsMouseOver" Value="True">
                                <Setter TargetName="border" Property="BorderBrush" Value="#C8D0E0"/>
                            </Trigger>
                        </ControlTemplate.Triggers>
                    </ControlTemplate>
                </Setter.Value>
            </Setter>
        </Style>

        <!-- Step indicator dot — hairline minimal su sfondo light -->
        <Style x:Key="StepDot" TargetType="Border">
            <Setter Property="Width" Value="26"/>
            <Setter Property="Height" Value="26"/>
            <Setter Property="CornerRadius" Value="13"/>
            <Setter Property="Background" Value="White"/>
            <Setter Property="BorderBrush" Value="#D6DCE8"/>
            <Setter Property="BorderThickness" Value="1"/>
        </Style>

        <!-- Card surface: hairline minimal panel light -->
        <Style x:Key="Card" TargetType="Border">
            <Setter Property="Background" Value="White"/>
            <Setter Property="BorderBrush" Value="#E4E8F0"/>
            <Setter Property="BorderThickness" Value="1"/>
            <Setter Property="CornerRadius" Value="6"/>
            <Setter Property="Padding" Value="20"/>
        </Style>
    </Window.Resources>

    <Grid>
        <Grid.ColumnDefinitions>
            <ColumnDefinition Width="280"/>
            <ColumnDefinition Width="*"/>
        </Grid.ColumnDefinitions>

        <!-- =============== SIDEBAR =============== -->
        <Grid Grid.Column="0" ClipToBounds="True" Background="#F4F6FB">
            <Grid.RowDefinitions>
                <RowDefinition Height="Auto"/>
                <RowDefinition Height="*"/>
                <RowDefinition Height="Auto"/>
            </Grid.RowDefinitions>

            <!-- Hairline divider verticale tra sidebar e content -->
            <Border Grid.Row="0" Grid.RowSpan="3" Width="1" Background="#E4E8F0" HorizontalAlignment="Right"/>

            <!-- Pattern decorativo: scatter network (richiama il logo). Punti magenta/blu del logo connessi da linee -->
            <Canvas Grid.Row="0" Grid.RowSpan="3" IsHitTestVisible="False" Opacity="0.5">
                <Line X1="50" Y1="220" X2="120" Y2="260" Stroke="#25D1397E" StrokeThickness="1"/>
                <Line X1="120" Y1="260" X2="180" Y2="230" Stroke="#25D1397E" StrokeThickness="1"/>
                <Line X1="120" Y1="260" X2="90" Y2="330" Stroke="#25D1397E" StrokeThickness="1"/>
                <Line X1="90" Y1="330" X2="180" Y2="350" Stroke="#25D1397E" StrokeThickness="1"/>
                <Line X1="180" Y1="350" X2="220" Y2="290" Stroke="#25D1397E" StrokeThickness="1"/>
                <Line X1="220" Y1="290" X2="180" Y2="230" Stroke="#25D1397E" StrokeThickness="1"/>
                <Line X1="90" Y1="330" X2="60" Y2="420" Stroke="#251A4FE0" StrokeThickness="1"/>
                <Line X1="60" Y1="420" X2="180" Y2="350" Stroke="#251A4FE0" StrokeThickness="1"/>

                <!-- Scatter points: cerchi piccoli con accent del logo -->
                <Ellipse Width="6" Height="6" Canvas.Left="47" Canvas.Top="217" Fill="#D1397E" Opacity="0.7"/>
                <Ellipse Width="6" Height="6" Canvas.Left="117" Canvas.Top="257" Fill="#D1397E" Opacity="0.7"/>
                <Ellipse Width="6" Height="6" Canvas.Left="177" Canvas.Top="227" Fill="#D1397E" Opacity="0.7"/>
                <Ellipse Width="6" Height="6" Canvas.Left="87" Canvas.Top="327" Fill="#1A4FE0" Opacity="0.7"/>
                <Ellipse Width="6" Height="6" Canvas.Left="177" Canvas.Top="347" Fill="#D1397E" Opacity="0.7"/>
                <Ellipse Width="6" Height="6" Canvas.Left="217" Canvas.Top="287" Fill="#D1397E" Opacity="0.7"/>
                <Ellipse Width="4" Height="4" Canvas.Left="58" Canvas.Top="418" Fill="#1A4FE0" Opacity="0.7"/>
            </Canvas>

            <!-- Logo PHASE in card minimal con hairline shadow soft -->
            <Border Grid.Row="0" Background="White" CornerRadius="10" Margin="24,32,24,24" Padding="18,16" BorderBrush="#E4E8F0" BorderThickness="1">
                <Border.Effect>
                    <DropShadowEffect Color="#FF1A1F2E" BlurRadius="12" ShadowDepth="2" Opacity="0.06"/>
                </Border.Effect>
                <Image x:Name="LogoImage" Stretch="Uniform" MaxHeight="78"/>
            </Border>

            <!-- Step indicator: hairline minimal, numerazione monospace -->
            <Grid Grid.Row="1" Margin="34,8,24,0">
                <!-- Linea verticale di connessione (hairline) -->
                <Border Background="#D6DCE8" Width="1" HorizontalAlignment="Left" Margin="13,18,0,12"/>
                <StackPanel>
                    <TextBlock Text="SETUP · PIPELINE" FontFamily="JetBrains Mono, Cascadia Code, Consolas" FontSize="9" FontWeight="SemiBold" Foreground="#8C95B8" Margin="0,0,0,18"/>

                    <Grid x:Name="StepRow1" Margin="0,5">
                        <Grid.ColumnDefinitions><ColumnDefinition Width="Auto"/><ColumnDefinition Width="*"/></Grid.ColumnDefinitions>
                        <Border x:Name="StepDot1" Style="{StaticResource StepDot}"><TextBlock x:Name="StepNum1" Text="01" FontFamily="JetBrains Mono, Cascadia Code, Consolas" FontSize="9" FontWeight="SemiBold" HorizontalAlignment="Center" VerticalAlignment="Center"/></Border>
                        <TextBlock x:Name="StepText1" Grid.Column="1" Text="Welcome" FontSize="13" VerticalAlignment="Center" Margin="18,0,0,0"/>
                    </Grid>
                    <Grid x:Name="StepRow2" Margin="0,5">
                        <Grid.ColumnDefinitions><ColumnDefinition Width="Auto"/><ColumnDefinition Width="*"/></Grid.ColumnDefinitions>
                        <Border x:Name="StepDot2" Style="{StaticResource StepDot}"><TextBlock x:Name="StepNum2" Text="02" FontFamily="JetBrains Mono, Cascadia Code, Consolas" FontSize="9" FontWeight="SemiBold" HorizontalAlignment="Center" VerticalAlignment="Center"/></Border>
                        <TextBlock x:Name="StepText2" Grid.Column="1" Text="MATLAB" FontSize="13" VerticalAlignment="Center" Margin="18,0,0,0"/>
                    </Grid>
                    <Grid x:Name="StepRow3" Margin="0,5">
                        <Grid.ColumnDefinitions><ColumnDefinition Width="Auto"/><ColumnDefinition Width="*"/></Grid.ColumnDefinitions>
                        <Border x:Name="StepDot3" Style="{StaticResource StepDot}"><TextBlock x:Name="StepNum3" Text="03" FontFamily="JetBrains Mono, Cascadia Code, Consolas" FontSize="9" FontWeight="SemiBold" HorizontalAlignment="Center" VerticalAlignment="Center"/></Border>
                        <TextBlock x:Name="StepText3" Grid.Column="1" Text="SNAP" FontSize="13" VerticalAlignment="Center" Margin="18,0,0,0"/>
                    </Grid>
                    <Grid x:Name="StepRow4" Margin="0,5">
                        <Grid.ColumnDefinitions><ColumnDefinition Width="Auto"/><ColumnDefinition Width="*"/></Grid.ColumnDefinitions>
                        <Border x:Name="StepDot4" Style="{StaticResource StepDot}"><TextBlock x:Name="StepNum4" Text="04" FontFamily="JetBrains Mono, Cascadia Code, Consolas" FontSize="9" FontWeight="SemiBold" HorizontalAlignment="Center" VerticalAlignment="Center"/></Border>
                        <TextBlock x:Name="StepText4" Grid.Column="1" Text="Python" FontSize="13" VerticalAlignment="Center" Margin="18,0,0,0"/>
                    </Grid>
                    <Grid x:Name="StepRow5" Margin="0,5">
                        <Grid.ColumnDefinitions><ColumnDefinition Width="Auto"/><ColumnDefinition Width="*"/></Grid.ColumnDefinitions>
                        <Border x:Name="StepDot5" Style="{StaticResource StepDot}"><TextBlock x:Name="StepNum5" Text="05" FontFamily="JetBrains Mono, Cascadia Code, Consolas" FontSize="9" FontWeight="SemiBold" HorizontalAlignment="Center" VerticalAlignment="Center"/></Border>
                        <TextBlock x:Name="StepText5" Grid.Column="1" Text="Destination" FontSize="13" VerticalAlignment="Center" Margin="18,0,0,0"/>
                    </Grid>
                    <Grid x:Name="StepRow6" Margin="0,5">
                        <Grid.ColumnDefinitions><ColumnDefinition Width="Auto"/><ColumnDefinition Width="*"/></Grid.ColumnDefinitions>
                        <Border x:Name="StepDot6" Style="{StaticResource StepDot}"><TextBlock x:Name="StepNum6" Text="06" FontFamily="JetBrains Mono, Cascadia Code, Consolas" FontSize="9" FontWeight="SemiBold" HorizontalAlignment="Center" VerticalAlignment="Center"/></Border>
                        <TextBlock x:Name="StepText6" Grid.Column="1" Text="Installation" FontSize="13" VerticalAlignment="Center" Margin="18,0,0,0"/>
                    </Grid>
                    <Grid x:Name="StepRow7" Margin="0,5">
                        <Grid.ColumnDefinitions><ColumnDefinition Width="Auto"/><ColumnDefinition Width="*"/></Grid.ColumnDefinitions>
                        <Border x:Name="StepDot7" Style="{StaticResource StepDot}"><TextBlock x:Name="StepNum7" Text="07" FontFamily="JetBrains Mono, Cascadia Code, Consolas" FontSize="9" FontWeight="SemiBold" HorizontalAlignment="Center" VerticalAlignment="Center"/></Border>
                        <TextBlock x:Name="StepText7" Grid.Column="1" Text="Finish" FontSize="13" VerticalAlignment="Center" Margin="18,0,0,0"/>
                    </Grid>
                </StackPanel>
            </Grid>

            <!-- Footer sidebar: tagline + build info, all monospace tone, light theme -->
            <StackPanel Grid.Row="2" Margin="34,0,24,28">
                <Border Background="White" CornerRadius="4" Padding="14,12" Margin="0,0,0,14" BorderBrush="#E4E8F0" BorderThickness="1">
                    <StackPanel>
                        <TextBlock FontSize="10" Foreground="#4A5168" Margin="0,0,0,4">
                            <Run Text="PERSISTENT  SCATTERER" FontWeight="Bold" Foreground="#0F1430"/>
                        </TextBlock>
                        <TextBlock Text="Highly Automated Suite for Environmental Monitoring" TextWrapping="Wrap" FontSize="10" Foreground="#8C95B8" LineHeight="14"/>
                    </StackPanel>
                </Border>
                <TextBlock FontFamily="JetBrains Mono, Cascadia Code, Consolas" FontSize="10" Foreground="#8C95B8">
                    <Run Text="build" Foreground="#B0B8C8"/>
                    <Run Text="  1.9.0  " Foreground="#1A4FE0" FontWeight="SemiBold"/>
                    <Run Text="·  pyccino/PHASE" Foreground="#8C95B8"/>
                </TextBlock>
            </StackPanel>
        </Grid>

        <!-- =============== CONTENT AREA =============== -->
        <Grid Grid.Column="1">
            <Grid.RowDefinitions>
                <RowDefinition Height="Auto"/>
                <RowDefinition Height="*"/>
                <RowDefinition Height="Auto"/>
            </Grid.RowDefinitions>

            <!-- Step label header — monospace, blu primary, breadcrumb tecnico -->
            <Border Grid.Row="0" Background="White" BorderBrush="#E4E8F0" BorderThickness="0,0,0,1" Padding="40,20">
                <Grid>
                    <Grid.ColumnDefinitions>
                        <ColumnDefinition Width="*"/>
                        <ColumnDefinition Width="Auto"/>
                    </Grid.ColumnDefinitions>
                    <TextBlock x:Name="StepLabel" FontFamily="JetBrains Mono, Cascadia Code, Consolas" FontSize="10" FontWeight="SemiBold" Foreground="#1A4FE0"/>
                    <TextBlock Grid.Column="1" FontFamily="JetBrains Mono, Cascadia Code, Consolas" FontSize="10" Foreground="#B0B8C8">
                        <Run Text="//" Foreground="#1A4FE0"/>
                        <Run Text=" phase-installer "/>
                        <Run Text="//" Foreground="#1A4FE0"/>
                        <Run Text=" pyccino"/>
                    </TextBlock>
                </Grid>
            </Border>

            <!-- Wizard pages -->
            <Grid Grid.Row="1" Margin="40,32,40,16">

            <!-- Page 1: Welcome -->
            <StackPanel x:Name="Page1_Welcome" Visibility="Visible">
                <TextBlock Text="PHASE Windows Installer" FontSize="28" FontWeight="Light" Margin="0,0,0,10"/>
                <TextBlock Text="End-to-end installation of the PHASE suite (Persistent scatterer Highly Automated Suite for Environmental monitoring) and all its dependencies on Windows."
                           TextWrapping="Wrap" FontSize="13" Foreground="#4A5168" Margin="0,0,0,24"/>
                <Border Style="{StaticResource Card}">
                    <StackPanel>
                        <TextBlock Text="What this wizard does" FontWeight="SemiBold" FontSize="13" Margin="0,0,0,14"/>
                        <TextBlock Text="1.  Detect MATLAB (must be installed and activated)" Margin="0,4"/>
                        <TextBlock Text="2.  Detect SNAP 13 (installs the bundled installer if missing)" Margin="0,4"/>
                        <TextBlock Text="3.  Install Python 3.11+ silently if missing" Margin="0,4"/>
                        <TextBlock Text="4.  Clone PHASE, StaMPS, TRAIN into the chosen folder" Margin="0,4"/>
                        <TextBlock Text="5.  Download native binaries, install GMT, configure all paths" Margin="0,4"/>
                    </StackPanel>
                </Border>
                <TextBlock Text="Estimated time: 15-30 minutes (depending on connection and SNAP installer)."
                           Margin="0,18,0,0" FontStyle="Italic" Foreground="#8C95B8" FontSize="12"/>
            </StackPanel>

            <!-- Page 2: MATLAB -->
            <StackPanel x:Name="Page2_Matlab" Visibility="Collapsed">
                <TextBlock Text="MATLAB" FontSize="28" FontWeight="Light" Margin="0,0,0,10"/>
                <TextBlock Text="PHASE requires MATLAB R2023a or newer. The installer cannot install MATLAB automatically because it is proprietary — it must already be installed and activated on this system."
                           TextWrapping="Wrap" FontSize="13" Foreground="#4A5168" Margin="0,0,0,20"/>

                <TextBlock x:Name="MatlabStatus" Text="" FontWeight="SemiBold" Margin="0,0,0,10"/>

                <TextBlock Text="Path to matlab.exe" Margin="0,0,0,6" FontSize="11" Foreground="#8C95B8" FontFamily="JetBrains Mono, Cascadia Code, Consolas"/>
                <Grid>
                    <Grid.ColumnDefinitions>
                        <ColumnDefinition Width="*"/>
                        <ColumnDefinition Width="Auto"/>
                    </Grid.ColumnDefinitions>
                    <TextBox x:Name="MatlabPathBox" Grid.Column="0"/>
                    <Button x:Name="MatlabBrowseBtn" Grid.Column="1" Content="Browse" Margin="8,0,0,0"/>
                </Grid>

                <TextBlock x:Name="MatlabHint" Text="" TextWrapping="Wrap" Margin="0,12,0,0" Foreground="#8C95B8" FontSize="12"/>

                <Border x:Name="MatlabDownloadHint" BorderBrush="#F0C492" BorderThickness="1" CornerRadius="4"
                        Padding="14,12" Background="#FFFBF1" Margin="0,16,0,0" Visibility="Collapsed">
                    <StackPanel>
                        <TextBlock Text="MATLAB not found on this system." FontWeight="SemiBold" Foreground="#A55B00"/>
                        <TextBlock Text="Download and install it from mathworks.com, then come back and provide the path manually."
                                   TextWrapping="Wrap" Margin="0,4,0,10" Foreground="#4A5168"/>
                        <Button x:Name="OpenMathworksBtn" Content="Open mathworks.com" HorizontalAlignment="Left"/>
                    </StackPanel>
                </Border>

                <Border x:Name="MatlabToolboxStatus" Style="{StaticResource Card}" Margin="0,16,0,0" Visibility="Collapsed">
                    <StackPanel>
                        <TextBlock x:Name="MatlabToolboxHeader" Text="MATLAB toolboxes" FontWeight="SemiBold" Margin="0,0,0,8"/>
                        <TextBlock x:Name="MatlabToolboxList" Text="" TextWrapping="Wrap" FontSize="12" FontFamily="JetBrains Mono, Cascadia Code, Consolas" LineHeight="20"/>
                        <TextBlock x:Name="MatlabToolboxHint" Text="" TextWrapping="Wrap" FontSize="11"
                                   Foreground="#8C95B8" Margin="0,10,0,0" Visibility="Collapsed"/>
                    </StackPanel>
                </Border>
            </StackPanel>

            <!-- Page 3: SNAP -->
            <StackPanel x:Name="Page3_Snap" Visibility="Collapsed">
                <TextBlock Text="SNAP" FontSize="28" FontWeight="Light" Margin="0,0,0,10"/>
                <TextBlock Text="ESA's Sentinel Application Platform (SNAP) is required for SLC preprocessing. Version 13.x is recommended."
                           TextWrapping="Wrap" FontSize="13" Foreground="#4A5168" Margin="0,0,0,20"/>

                <TextBlock x:Name="SnapStatus" Text="" FontWeight="SemiBold" Margin="0,0,0,10"/>

                <TextBlock Text="Path to gpt.exe" Margin="0,0,0,6" FontSize="11" Foreground="#8C95B8" FontFamily="JetBrains Mono, Cascadia Code, Consolas"/>
                <Grid>
                    <Grid.ColumnDefinitions>
                        <ColumnDefinition Width="*"/>
                        <ColumnDefinition Width="Auto"/>
                    </Grid.ColumnDefinitions>
                    <TextBox x:Name="SnapPathBox" Grid.Column="0"/>
                    <Button x:Name="SnapBrowseBtn" Grid.Column="1" Content="Browse" Margin="8,0,0,0"/>
                </Grid>

                <Border x:Name="SnapInstallHint" BorderBrush="#F0C492" BorderThickness="1" CornerRadius="4"
                        Padding="14,12" Background="#FFFBF1" Margin="0,16,0,0" Visibility="Collapsed">
                    <StackPanel>
                        <TextBlock Text="SNAP not found on this system." FontWeight="SemiBold" Foreground="#A55B00"/>
                        <TextBlock x:Name="SnapInstallText" Text="" TextWrapping="Wrap" Margin="0,4,0,10" Foreground="#4A5168"/>
                        <Button x:Name="InstallSnapBtn" Content="Install SNAP now" HorizontalAlignment="Left"/>
                        <ProgressBar x:Name="SnapProgress" Height="6" Margin="0,14,0,4" Visibility="Collapsed" Maximum="100" Background="#E4E8F0" Foreground="#1A4FE0" BorderThickness="0"/>
                        <TextBlock x:Name="SnapProgressText" Text="" FontSize="11" Foreground="#4A5168" FontFamily="JetBrains Mono, Cascadia Code, Consolas" Visibility="Collapsed"/>
                    </StackPanel>
                </Border>

                <TextBlock x:Name="SnapHint" Text="" TextWrapping="Wrap" Margin="0,12,0,0" Foreground="#8C95B8" FontSize="12"/>
            </StackPanel>

            <!-- Page 4: Python -->
            <StackPanel x:Name="Page4_Python" Visibility="Collapsed">
                <TextBlock Text="Python" FontSize="28" FontWeight="Light" Margin="0,0,0,10"/>
                <TextBlock Text="PHASE requires Python 3.11 or newer with the openpyxl, requests, asf_search and shapely libraries. The installer can download and install everything automatically."
                           TextWrapping="Wrap" FontSize="13" Foreground="#4A5168" Margin="0,0,0,20"/>

                <TextBlock x:Name="PythonStatus" Text="" FontWeight="SemiBold" Margin="0,0,0,12"/>
                <TextBlock x:Name="PythonPath" Text="" FontFamily="JetBrains Mono, Cascadia Code, Consolas" FontSize="12" TextWrapping="Wrap" Margin="0,0,0,12" Foreground="#4A5168"/>

                <Border x:Name="PythonInstallPanel" BorderBrush="#F0C492" BorderThickness="1" CornerRadius="4"
                        Padding="14,12" Background="#FFFBF1" Margin="0,8,0,0" Visibility="Collapsed">
                    <StackPanel>
                        <TextBlock Text="Python 3.11+ not detected (or only the Microsoft Store stub was found)." FontWeight="SemiBold" Foreground="#A55B00"/>
                        <TextBlock Text="It will be downloaded from python.org and installed silently for the current user (~28 MB)."
                                   TextWrapping="Wrap" Margin="0,4,0,10" Foreground="#4A5168"/>
                        <Button x:Name="InstallPythonBtn" Content="Install Python 3.11.9 now" HorizontalAlignment="Left"/>
                        <ProgressBar x:Name="PythonProgress" Height="6" Margin="0,14,0,4" Visibility="Collapsed" Maximum="100" Background="#E4E8F0" Foreground="#1A4FE0" BorderThickness="0"/>
                        <TextBlock x:Name="PythonProgressText" Text="" FontSize="11" Foreground="#4A5168" FontFamily="JetBrains Mono, Cascadia Code, Consolas" Visibility="Collapsed"/>
                    </StackPanel>
                </Border>
            </StackPanel>

            <!-- Page 5: Destination folder -->
            <StackPanel x:Name="Page5_Dest" Visibility="Collapsed">
                <TextBlock Text="Destination folder" FontSize="28" FontWeight="Light" Margin="0,0,0,10"/>
                <TextBlock Text="Choose where PHASE will be installed. Three subfolders will be created: PHASE\, StaMPS\, TRAIN\."
                           TextWrapping="Wrap" FontSize="13" Foreground="#4A5168" Margin="0,0,0,20"/>

                <TextBlock Text="Folder" Margin="0,0,0,6" FontSize="11" Foreground="#8C95B8" FontFamily="JetBrains Mono, Cascadia Code, Consolas"/>
                <Grid>
                    <Grid.ColumnDefinitions>
                        <ColumnDefinition Width="*"/>
                        <ColumnDefinition Width="Auto"/>
                    </Grid.ColumnDefinitions>
                    <TextBox x:Name="DestPathBox" Grid.Column="0"/>
                    <Button x:Name="DestBrowseBtn" Grid.Column="1" Content="Browse" Margin="8,0,0,0"/>
                </Grid>

                <TextBlock x:Name="DestHint" Text="" TextWrapping="Wrap" Margin="0,12,0,0" Foreground="#8C95B8" FontSize="12"/>

                <Border Style="{StaticResource Card}" Margin="0,18,0,0">
                    <StackPanel>
                        <TextBlock Text="Recommendations" FontWeight="SemiBold" Margin="0,0,0,8"/>
                        <TextBlock Text="·  Short paths (near drive root) avoid the MAX_PATH=260 limit on StaMPS PATCH_N/ directories." Margin="0,3" Foreground="#4A5168" TextWrapping="Wrap"/>
                        <TextBlock Text="·  Avoid folders inside OneDrive (intermittent file locks during long runs)." Margin="0,3" Foreground="#4A5168" TextWrapping="Wrap"/>
                        <TextBlock Text="·  ASCII characters only (no accents, spaces are fine)." Margin="0,3" Foreground="#4A5168" TextWrapping="Wrap"/>
                    </StackPanel>
                </Border>
            </StackPanel>

            <!-- Page 6: Installation — task list view (replaces console log) -->
            <StackPanel x:Name="Page6_Setup" Visibility="Collapsed">
                <Grid Margin="0,0,0,14">
                    <Grid.ColumnDefinitions>
                        <ColumnDefinition Width="*"/>
                        <ColumnDefinition Width="Auto"/>
                    </Grid.ColumnDefinitions>
                    <StackPanel Grid.Column="0">
                        <TextBlock Text="Installation" FontSize="28" FontWeight="Light" Margin="0,0,0,4"/>
                        <TextBlock x:Name="SetupSubtitle" Text="Ready to install. Click Start to begin."
                                   TextWrapping="Wrap" FontSize="13" Foreground="#4A5168"/>
                    </StackPanel>
                    <Button x:Name="StartSetupBtn" Grid.Column="1" Content="Start installation"
                            Style="{StaticResource PrimaryButton}" VerticalAlignment="Center"/>
                </Grid>

                <!-- Task list (replaces the old dark console log) -->
                <Border Style="{StaticResource Card}" Padding="22,18">
                    <StackPanel>
                        <Grid Margin="0,0,0,14">
                            <Grid.ColumnDefinitions>
                                <ColumnDefinition Width="*"/>
                                <ColumnDefinition Width="Auto"/>
                            </Grid.ColumnDefinitions>
                            <TextBlock Text="Pipeline" FontFamily="JetBrains Mono, Cascadia Code, Consolas" FontSize="10" FontWeight="SemiBold" Foreground="#1A4FE0"/>
                            <TextBlock Grid.Column="1" x:Name="TaskCounter" Text="0 / 10" FontFamily="JetBrains Mono, Cascadia Code, Consolas" FontSize="10" Foreground="#8C95B8"/>
                        </Grid>

                        <ScrollViewer x:Name="TaskListScroll" Height="300" VerticalScrollBarVisibility="Auto" HorizontalScrollBarVisibility="Disabled">
                            <StackPanel x:Name="TaskList"/>
                        </ScrollViewer>
                    </StackPanel>
                </Border>

                <!-- Bottom progress bar — slim, like loaders in modern devtools -->
                <Grid Margin="0,14,0,0">
                    <Grid.ColumnDefinitions>
                        <ColumnDefinition Width="*"/>
                        <ColumnDefinition Width="Auto"/>
                    </Grid.ColumnDefinitions>
                    <ProgressBar x:Name="SetupProgress" Height="4" Maximum="100"
                                 Background="#E4E8F0" Foreground="#1A4FE0" BorderThickness="0" VerticalAlignment="Center"/>
                    <TextBlock Grid.Column="1" x:Name="SetupProgressText" Text="" Margin="14,0,0,0"
                               Foreground="#4A5168" FontSize="11" FontFamily="JetBrains Mono, Cascadia Code, Consolas" VerticalAlignment="Center"/>
                </Grid>
            </StackPanel>

            <!-- Page 7: Finish -->
            <StackPanel x:Name="Page7_Finish" Visibility="Collapsed">
                <TextBlock Text="Installation complete" FontSize="28" FontWeight="Light" Foreground="#2DBA6E" Margin="0,0,0,10"/>
                <TextBlock x:Name="FinishSubtitle" Text="PHASE is ready. Open the folder and double-click one of the three .mlapp files."
                           TextWrapping="Wrap" FontSize="13" Foreground="#4A5168" Margin="0,0,0,22"/>

                <Border Style="{StaticResource Card}">
                    <StackPanel>
                        <TextBlock Text="PHASE folder" FontFamily="JetBrains Mono, Cascadia Code, Consolas" FontSize="10" FontWeight="SemiBold" Foreground="#1A4FE0" Margin="0,0,0,4"/>
                        <TextBlock x:Name="FinishPath" Text="" FontFamily="JetBrains Mono, Cascadia Code, Consolas" FontSize="12" Margin="0,0,0,18" Foreground="#0F1430"/>
                        <TextBlock Text="Available MATLAB apps" FontFamily="JetBrains Mono, Cascadia Code, Consolas" FontSize="10" FontWeight="SemiBold" Foreground="#1A4FE0" Margin="0,0,0,8"/>
                        <TextBlock Text="·  PHASE_Preprocessing.mlapp  —  module 1 (SNAP preprocessing)" Margin="0,3" Foreground="#4A5168"/>
                        <TextBlock Text="·  PHASE_Preprocessing\PHASE_StaMPS.mlapp  —  module 2 (PSI StaMPS)" Margin="0,3" Foreground="#4A5168"/>
                        <TextBlock Text="·  PHASE_model.mlapp  —  module 3 (geospatial analysis)" Margin="0,3" Foreground="#4A5168"/>
                    </StackPanel>
                </Border>

                <StackPanel Orientation="Horizontal" Margin="0,22,0,0">
                    <Button x:Name="OpenFolderBtn" Content="Open PHASE folder"/>
                    <Button x:Name="OpenLogBtn" Content="Open install log"/>
                </StackPanel>
            </StackPanel>

        </Grid>

            <!-- Footer with Back/Next/Cancel — light, hairline divider -->
            <Border Grid.Row="2" Background="#FCFCFE" BorderBrush="#E4E8F0" BorderThickness="0,1,0,0" Padding="40,18">
                <StackPanel Orientation="Horizontal" HorizontalAlignment="Right">
                    <Button x:Name="BackBtn" Content="Back" Style="{StaticResource SecondaryButton}"/>
                    <Button x:Name="NextBtn" Content="Next" Style="{StaticResource PrimaryButton}"/>
                    <Button x:Name="CancelBtn" Content="Cancel" Style="{StaticResource SecondaryButton}"/>
                </StackPanel>
            </Border>
        </Grid>
    </Grid>
</Window>
'@

# -----------------------------------------------------------------------------
# Parse XAML and wire up event handlers
# -----------------------------------------------------------------------------

$reader = New-Object System.Xml.XmlNodeReader $xaml
$window = [Windows.Markup.XamlReader]::Load($reader)

# Helper: find named element
function Get-Element { param([string]$Name) $window.FindName($Name) }

# Carica il logo PHASE nella sidebar. Lo cerchiamo accanto all'.exe (file
# bundled in installer/PHASE_logo.png) e graceful-degrade se assente.
function Set-PhaseLogo {
    $candidates = @(
        (Join-Path $Script:ScriptDir 'PHASE_logo.png'),
        (Join-Path $Script:ScriptDir '..\PHASE_logo.png'),
        'F:\phase\PHASE_logo.png'
    )
    foreach ($p in $candidates) {
        if ($p -and (Test-Path $p)) {
            try {
                $bmp = New-Object System.Windows.Media.Imaging.BitmapImage
                $bmp.BeginInit()
                $bmp.UriSource = New-Object System.Uri ($p, [System.UriKind]::Absolute)
                $bmp.CacheOption = 'OnLoad'
                $bmp.EndInit()
                $bmp.Freeze()
                (Get-Element 'LogoImage').Source = $bmp
                return
            } catch {
                continue
            }
        }
    }
    # Fallback: testo se nessun logo trovato (mantengo Image vuota, il
    # Border bianco resta ma con un placeholder testuale)
    (Get-Element 'LogoImage').Source = $null
}
Set-PhaseLogo

# Element references
$pages = @{
    1 = Get-Element 'Page1_Welcome'
    2 = Get-Element 'Page2_Matlab'
    3 = Get-Element 'Page3_Snap'
    4 = Get-Element 'Page4_Python'
    5 = Get-Element 'Page5_Dest'
    6 = Get-Element 'Page6_Setup'
    7 = Get-Element 'Page7_Finish'
}
$labels = @{
    1 = 'WELCOME'
    2 = 'STEP 2  ·  MATLAB'
    3 = 'STEP 3  ·  SNAP'
    4 = 'STEP 4  ·  PYTHON'
    5 = 'STEP 5  ·  DESTINATION'
    6 = 'STEP 6  ·  INSTALLATION'
    7 = 'STEP 7  ·  FINISH'
}

$Script:CurrentPage = 1
$Script:SetupLogPath = Join-Path $env:TEMP "phase-installer-$(Get-Date -Format 'yyyyMMdd-HHmmss').log"

# Aggiorna lo step indicator nella sidebar (light theme).
# done    = blu solido con check bianco
# current = blu solido con numero bianco + soft glow
# pending = bianco con hairline grigio, numero muted
function Update-StepIndicator {
    param([int]$Current)
    $blue        = [System.Windows.Media.ColorConverter]::ConvertFromString('#FF1A4FE0')
    $blueLight   = [System.Windows.Media.ColorConverter]::ConvertFromString('#FFE8EEFF')
    $white       = [System.Windows.Media.ColorConverter]::ConvertFromString('#FFFFFFFF')
    $borderLine  = [System.Windows.Media.ColorConverter]::ConvertFromString('#FFD6DCE8')
    $ink         = [System.Windows.Media.ColorConverter]::ConvertFromString('#FF0F1430')
    $textMuted   = [System.Windows.Media.ColorConverter]::ConvertFromString('#FF8C95B8')

    # Soft glow per il dot corrente
    $glowEffect = New-Object System.Windows.Media.Effects.DropShadowEffect
    $glowEffect.Color = [System.Windows.Media.ColorConverter]::ConvertFromString('#FF1A4FE0')
    $glowEffect.BlurRadius = 14
    $glowEffect.ShadowDepth = 0
    $glowEffect.Opacity = 0.35

    $blueBrush  = New-Object System.Windows.Media.SolidColorBrush $blue
    $whiteBrush = New-Object System.Windows.Media.SolidColorBrush $white
    $inkBrush   = New-Object System.Windows.Media.SolidColorBrush $ink
    $mutedBrush = New-Object System.Windows.Media.SolidColorBrush $textMuted
    $borderBrush = New-Object System.Windows.Media.SolidColorBrush $borderLine

    for ($i = 1; $i -le 7; $i++) {
        $dot = Get-Element "StepDot$i"
        $num = Get-Element "StepNum$i"
        $txt = Get-Element "StepText$i"
        if (-not $dot) { continue }
        if ($i -lt $Current) {
            # Done: dot blu pieno con check bianco
            $dot.Background = $blueBrush
            $dot.BorderBrush = $blueBrush
            $dot.Effect = $null
            $num.Text = [char]0x2713
            $num.Foreground = $whiteBrush
            $num.FontWeight = 'Bold'
            $txt.Foreground = $inkBrush
            $txt.Opacity = 0.85
            $txt.FontWeight = 'Normal'
        } elseif ($i -eq $Current) {
            # Current: dot blu pieno con numero bianco + glow
            $dot.Background = $blueBrush
            $dot.BorderBrush = $blueBrush
            $dot.Effect = $glowEffect
            $num.Text = ("{0:D2}" -f $i)
            $num.Foreground = $whiteBrush
            $num.FontWeight = 'Bold'
            $txt.Foreground = $inkBrush
            $txt.Opacity = 1.0
            $txt.FontWeight = 'Bold'
        } else {
            # Pending: dot bianco con hairline, numero muted
            $dot.Background = $whiteBrush
            $dot.BorderBrush = $borderBrush
            $dot.Effect = $null
            $num.Text = ("{0:D2}" -f $i)
            $num.Foreground = $mutedBrush
            $num.FontWeight = 'SemiBold'
            $txt.Foreground = $mutedBrush
            $txt.Opacity = 1.0
            $txt.FontWeight = 'Normal'
        }
    }
}

function Show-Page {
    param([int]$N)
    foreach ($k in $pages.Keys) {
        $pages[$k].Visibility = if ($k -eq $N) { 'Visible' } else { 'Collapsed' }
    }
    (Get-Element 'StepLabel').Text = $labels[$N]
    Update-StepIndicator -Current $N
    (Get-Element 'BackBtn').IsEnabled = ($N -gt 1 -and $N -lt 7)
    (Get-Element 'NextBtn').IsEnabled = ($N -lt 6)
    (Get-Element 'NextBtn').Content = if ($N -eq 5) { 'Continue' } else { 'Next' }
    if ($N -eq 7) {
        (Get-Element 'NextBtn').Visibility = 'Collapsed'
        (Get-Element 'CancelBtn').Content = 'Close'
    }
    $Script:CurrentPage = $N

    # Page-specific initialization
    switch ($N) {
        2 { Initialize-MatlabPage }
        3 { Initialize-SnapPage }
        4 { Initialize-PythonPage }
        5 { Initialize-DestPage }
        6 { Initialize-SetupPage }
        7 { Initialize-FinishPage }
    }
}

function Initialize-MatlabPage {
    $found = Find-Matlab
    if ($found) {
        (Get-Element 'MatlabStatus').Text = "[OK] MATLAB detected automatically."
        (Get-Element 'MatlabStatus').Foreground = '#FF008000'
        (Get-Element 'MatlabPathBox').Text = $found
        (Get-Element 'MatlabDownloadHint').Visibility = 'Collapsed'
        (Get-Element 'MatlabHint').Text = "Path detected automatically. You can override it to target a different MATLAB install."
        Update-MatlabToolboxStatus -MatlabExe $found
    } else {
        (Get-Element 'MatlabStatus').Text = "[X] MATLAB not found on system."
        (Get-Element 'MatlabStatus').Foreground = '#FFA04000'
        (Get-Element 'MatlabPathBox').Text = ''
        (Get-Element 'MatlabDownloadHint').Visibility = 'Visible'
        (Get-Element 'MatlabHint').Text = "Inserisci il path completo a matlab.exe (es. C:\Program Files\MATLAB\R2025a\bin\matlab.exe)."
        (Get-Element 'MatlabToolboxStatus').Visibility = 'Collapsed'
    }
    Update-MatlabValidation
}

# Aggiorna il pannello "Toolbox MATLAB" della Page 2 con detection via
# filesystem (<MATLABROOT>\toolbox\<name>\). Non-bloccante: mostra solo
# lo stato, l'utente prosegue comunque.
function Update-MatlabToolboxStatus {
    param([string]$MatlabExe)
    if (-not $MatlabExe -or -not (Test-Path $MatlabExe)) {
        (Get-Element 'MatlabToolboxStatus').Visibility = 'Collapsed'
        return
    }
    try {
        $toolboxes = Find-MatlabToolboxes -MatlabExe $MatlabExe
    } catch {
        (Get-Element 'MatlabToolboxStatus').Visibility = 'Collapsed'
        return
    }

    # Costruisco la lista visibile + flag se qualcosa di richiesto manca
    $lines = @()
    $missingRequired = @()
    foreach ($name in $toolboxes.Keys) {
        $installed = $toolboxes[$name]
        $isRequired = $Script:RequiredToolboxes -contains $name
        $marker = if ($installed) { '[OK]' } else { if ($isRequired) { '[X]' } else { '[ ]' } }
        $tag = if ($isRequired) { ' (richiesta)' } else { ' (opzionale)' }
        $lines += "$marker $name$tag"
        if (-not $installed -and $isRequired) { $missingRequired += $name }
    }

    (Get-Element 'MatlabToolboxList').Text = ($lines -join "`n")

    if ($missingRequired.Count -eq 0) {
        (Get-Element 'MatlabToolboxHeader').Text = "Toolbox MATLAB - tutte le richieste presenti"
        (Get-Element 'MatlabToolboxHeader').Foreground = '#FF008000'
        (Get-Element 'MatlabToolboxHint').Visibility = 'Collapsed'
    } else {
        (Get-Element 'MatlabToolboxHeader').Text = "Toolbox MATLAB - $($missingRequired.Count) richiesta/e mancanti"
        (Get-Element 'MatlabToolboxHeader').Foreground = '#FFA04000'
        (Get-Element 'MatlabToolboxHint').Visibility = 'Visible'
        (Get-Element 'MatlabToolboxHint').Text = "To install missing toolboxes: open MATLAB -> Home -> Add-Ons -> Get Add-Ons -> search by name -> Install. You are already signed in to MATLAB, no extra credentials required. You can continue the wizard now and add the toolboxes later."
    }
    (Get-Element 'MatlabToolboxStatus').Visibility = 'Visible'
}

function Update-MatlabValidation {
    $path = (Get-Element 'MatlabPathBox').Text
    $valid = ($path -and (Test-Path $path) -and ($path -like '*matlab.exe'))
    (Get-Element 'NextBtn').IsEnabled = $valid
    if ($valid) {
        $Script:State.MatlabExe = $path
        Update-MatlabToolboxStatus -MatlabExe $path
    } else {
        (Get-Element 'MatlabToolboxStatus').Visibility = 'Collapsed'
    }
}

function Initialize-SnapPage {
    $found = Find-Snap
    if ($found) {
        (Get-Element 'SnapStatus').Text = "[OK] SNAP detected automatically."
        (Get-Element 'SnapStatus').Foreground = '#FF008000'
        (Get-Element 'SnapPathBox').Text = $found
        (Get-Element 'SnapInstallHint').Visibility = 'Collapsed'
        (Get-Element 'SnapHint').Text = "Path detected automatically. Override it to point at a different install."
    } else {
        (Get-Element 'SnapStatus').Text = "[X] SNAP not found on system."
        (Get-Element 'SnapStatus').Foreground = '#FFA04000'
        (Get-Element 'SnapPathBox').Text = ''
        $hint = if (Test-Path $Script:BundledSnapPath) {
            "SNAP 13.0.0 installer is bundled with this distribution (~1 GB). Clicking 'Install SNAP now' runs it silently: a single UAC prompt is required, ~3-5 minutes total. Silent failure falls back to the interactive wizard automatically."
        } else {
            "SNAP 13.0.0 will be downloaded from ESA (~1 GB, 2-5 min on a fast connection) and installed silently. Only a single UAC prompt is required during the install step."
        }
        (Get-Element 'SnapInstallText').Text = $hint
        (Get-Element 'SnapInstallHint').Visibility = 'Visible'
        # Always enabled — we can either use the bundle or download on the fly
        (Get-Element 'InstallSnapBtn').IsEnabled = $true
    }
    Update-SnapValidation
}

function Update-SnapValidation {
    $path = (Get-Element 'SnapPathBox').Text
    $valid = ($path -and (Test-Path $path) -and ($path -like '*gpt.exe'))
    (Get-Element 'NextBtn').IsEnabled = $valid
    if ($valid) { $Script:State.SnapGpt = $path }
}

function Initialize-PythonPage {
    $found = Find-Python
    if ($found) {
        (Get-Element 'PythonStatus').Text = "[OK] Python $($found.Version) detected."
        (Get-Element 'PythonStatus').Foreground = '#FF008000'
        (Get-Element 'PythonPath').Text = $found.Exe
        (Get-Element 'PythonInstallPanel').Visibility = 'Collapsed'
        $Script:State.PythonExe = $found.Exe
        $Script:State.PythonVersion = $found.Version
        (Get-Element 'NextBtn').IsEnabled = $true
    } else {
        (Get-Element 'PythonStatus').Text = "[X] Python 3.11+ not detected."
        (Get-Element 'PythonStatus').Foreground = '#FFA04000'
        (Get-Element 'PythonPath').Text = ''
        (Get-Element 'PythonInstallPanel').Visibility = 'Visible'
        (Get-Element 'NextBtn').IsEnabled = $false
    }
}

function Initialize-DestPage {
    if (-not (Get-Element 'DestPathBox').Text) {
        (Get-Element 'DestPathBox').Text = $Script:State.InstallDir
    }
    Update-DestValidation
}

function Update-DestValidation {
    $path = (Get-Element 'DestPathBox').Text
    $valid = $false
    $hint = ''
    if (-not $path) {
        $hint = 'Specify a destination folder.'
    } elseif ($path -match '[^\x00-\x7F]') {
        $hint = 'Non-ASCII characters detected. Use a folder with ASCII characters only.'
    } elseif ($path -like '*OneDrive*') {
        $hint = '[!] Path under OneDrive: risk of file locks during long runs. Recommended to change.'
        $valid = $true   # warning, non blocco
    } else {
        $parent = Split-Path -Parent $path
        if (-not $parent -or (Test-Path $parent)) {
            $valid = $true
            $hint = if (Test-Path $path) { "A 'PHASE' folder will be created inside this folder (reused if it already exists)." } else { "This folder will be created, with a 'PHASE' folder inside it." }
        } else {
            $hint = "The parent folder $parent does not exist."
        }
    }
    (Get-Element 'DestHint').Text = $hint
    (Get-Element 'NextBtn').IsEnabled = $valid
    if ($valid) { $Script:State.InstallDir = $path }
}

# Definition of the 10 pipeline tasks shown in the Step 6 task list.
# These are visual checkpoints — they roughly mirror Invoke-FullSetup
# stages but their granularity is presentational, not 1:1 with code.
$Script:PipelineTasks = @(
    @{ Key = 'git';         Label = 'Verify Git installation' }
    @{ Key = 'clone-phase'; Label = 'Clone PHASE repository' }
    @{ Key = 'clone-stamps'; Label = 'Clone StaMPS repository' }
    @{ Key = 'clone-train'; Label = 'Clone TRAIN repository' }
    @{ Key = 'stamps-bin';  Label = 'Download native StaMPS binaries (precompiled)' }
    @{ Key = 'gmt';         Label = 'Install GMT (portable)' }
    @{ Key = 'env';         Label = 'Configure environment variables' }
    @{ Key = 'matlab';      Label = 'MATLAB savepath + precompile .mat files' }
    @{ Key = 'patch';       Label = 'Patch .mlapp files for auto-load' }
)

# Tracks the start time of each running task so we can report elapsed time
$Script:TaskTimers = @{}

function Initialize-SetupPage {
    Reset-TaskList
    (Get-Element 'SetupProgress').Value = 0
    (Get-Element 'SetupProgressText').Text = 'idle'
    (Get-Element 'BackBtn').IsEnabled = $false
    (Get-Element 'NextBtn').IsEnabled = $false
    (Get-Element 'StartSetupBtn').IsEnabled = $true
    (Get-Element 'TaskCounter').Text = ("0 / {0}" -f $Script:PipelineTasks.Count)
}

function Initialize-FinishPage {
    (Get-Element 'FinishPath').Text = (Join-Path $Script:State.InstallDir 'PHASE')
    (Get-Element 'BackBtn').IsEnabled = $false
    (Get-Element 'CancelBtn').IsEnabled = $true
    (Get-Element 'CancelBtn').Content = 'Close'
}

# Rebuilds the task list panel from scratch (used on page enter).
# Each task is rendered as a 3-column row: status icon | label+detail | duration.
function Reset-TaskList {
    $panel = Get-Element 'TaskList'
    $panel.Children.Clear()
    $Script:TaskNodes = @{}
    $Script:TaskTimers = @{}

    foreach ($task in $Script:PipelineTasks) {
        $row = New-Object System.Windows.Controls.Grid
        $row.Margin = '0,5'
        $col1 = New-Object System.Windows.Controls.ColumnDefinition; $col1.Width = 'Auto'
        $col2 = New-Object System.Windows.Controls.ColumnDefinition; $col2.Width = '*'
        $col3 = New-Object System.Windows.Controls.ColumnDefinition; $col3.Width = 'Auto'
        [void]$row.ColumnDefinitions.Add($col1)
        [void]$row.ColumnDefinitions.Add($col2)
        [void]$row.ColumnDefinitions.Add($col3)

        # Status pill: a 18 px circle Border with a TextBlock glyph inside
        $dot = New-Object System.Windows.Controls.Border
        $dot.Width = 18; $dot.Height = 18
        $dot.CornerRadius = New-Object System.Windows.CornerRadius 9
        $dot.Background = New-Object System.Windows.Media.SolidColorBrush ([System.Windows.Media.ColorConverter]::ConvertFromString('#FFFFFFFF'))
        $dot.BorderBrush = New-Object System.Windows.Media.SolidColorBrush ([System.Windows.Media.ColorConverter]::ConvertFromString('#FFD6DCE8'))
        $dot.BorderThickness = New-Object System.Windows.Thickness 1
        $dot.VerticalAlignment = 'Center'
        [System.Windows.Controls.Grid]::SetColumn($dot, 0)

        $glyph = New-Object System.Windows.Controls.TextBlock
        $glyph.Text = ''
        $glyph.HorizontalAlignment = 'Center'
        $glyph.VerticalAlignment = 'Center'
        $glyph.FontSize = 10
        $glyph.FontWeight = 'Bold'
        $glyph.Foreground = New-Object System.Windows.Media.SolidColorBrush ([System.Windows.Media.ColorConverter]::ConvertFromString('#FF8C95B8'))
        $dot.Child = $glyph

        # Center column: label + detail (stacked)
        $textCol = New-Object System.Windows.Controls.StackPanel
        $textCol.Margin = '14,0,12,0'
        $textCol.VerticalAlignment = 'Center'
        [System.Windows.Controls.Grid]::SetColumn($textCol, 1)

        $label = New-Object System.Windows.Controls.TextBlock
        $label.Text = $task.Label
        $label.FontSize = 13
        $label.Foreground = New-Object System.Windows.Media.SolidColorBrush ([System.Windows.Media.ColorConverter]::ConvertFromString('#FF8C95B8'))
        [void]$textCol.Children.Add($label)

        $detail = New-Object System.Windows.Controls.TextBlock
        $detail.Text = ''
        $detail.FontFamily = New-Object System.Windows.Media.FontFamily 'JetBrains Mono, Cascadia Code, Consolas'
        $detail.FontSize = 10
        $detail.Foreground = New-Object System.Windows.Media.SolidColorBrush ([System.Windows.Media.ColorConverter]::ConvertFromString('#FF8C95B8'))
        $detail.Margin = '0,2,0,0'
        $detail.Visibility = 'Collapsed'
        $detail.TextTrimming = 'CharacterEllipsis'
        [void]$textCol.Children.Add($detail)

        # Right column: duration display (mono)
        $duration = New-Object System.Windows.Controls.TextBlock
        $duration.Text = ''
        $duration.FontFamily = New-Object System.Windows.Media.FontFamily 'JetBrains Mono, Cascadia Code, Consolas'
        $duration.FontSize = 10
        $duration.Foreground = New-Object System.Windows.Media.SolidColorBrush ([System.Windows.Media.ColorConverter]::ConvertFromString('#FFB0B8C8'))
        $duration.VerticalAlignment = 'Center'
        [System.Windows.Controls.Grid]::SetColumn($duration, 2)

        [void]$row.Children.Add($dot)
        [void]$row.Children.Add($textCol)
        [void]$row.Children.Add($duration)

        [void]$panel.Children.Add($row)

        $Script:TaskNodes[$task.Key] = @{
            Row = $row; Dot = $dot; Glyph = $glyph;
            Label = $label; Detail = $detail; Duration = $duration
        }
    }
}

# Update a single task's visual state. Status: 'pending' | 'running' | 'done' | 'skip' | 'error'.
function Update-Task {
    param(
        [Parameter(Mandatory)] [string]$Key,
        [Parameter(Mandatory)] [ValidateSet('pending','running','done','skip','error')] [string]$Status,
        [string]$Detail
    )
    $node = $Script:TaskNodes[$Key]
    if (-not $node) { return }

    $ink   = [System.Windows.Media.ColorConverter]::ConvertFromString('#FF0F1430')
    $blue  = [System.Windows.Media.ColorConverter]::ConvertFromString('#FF1A4FE0')
    $green = [System.Windows.Media.ColorConverter]::ConvertFromString('#FF2DBA6E')
    $red   = [System.Windows.Media.ColorConverter]::ConvertFromString('#FFE03B5C')
    $muted = [System.Windows.Media.ColorConverter]::ConvertFromString('#FF8C95B8')
    $line  = [System.Windows.Media.ColorConverter]::ConvertFromString('#FFD6DCE8')
    $white = [System.Windows.Media.ColorConverter]::ConvertFromString('#FFFFFFFF')

    $glowBlue = New-Object System.Windows.Media.Effects.DropShadowEffect
    $glowBlue.Color = $blue; $glowBlue.BlurRadius = 10; $glowBlue.ShadowDepth = 0; $glowBlue.Opacity = 0.45

    switch ($Status) {
        'pending' {
            $node.Dot.Background = New-Object System.Windows.Media.SolidColorBrush $white
            $node.Dot.BorderBrush = New-Object System.Windows.Media.SolidColorBrush $line
            $node.Dot.Effect = $null
            $node.Glyph.Text = ''
            $node.Label.Foreground = New-Object System.Windows.Media.SolidColorBrush $muted
            $node.Label.FontWeight = 'Normal'
            $node.Detail.Visibility = 'Collapsed'
            $node.Duration.Text = ''
        }
        'running' {
            $node.Dot.Background = New-Object System.Windows.Media.SolidColorBrush $blue
            $node.Dot.BorderBrush = New-Object System.Windows.Media.SolidColorBrush $blue
            $node.Dot.Effect = $glowBlue
            $node.Glyph.Text = [char]0x25CF   # filled circle (will pulse via storyboard)
            $node.Glyph.Foreground = New-Object System.Windows.Media.SolidColorBrush $white
            $node.Label.Foreground = New-Object System.Windows.Media.SolidColorBrush $ink
            $node.Label.FontWeight = 'SemiBold'
            $Script:TaskTimers[$Key] = Get-Date
            # Pulse animation on the dot
            $anim = New-Object System.Windows.Media.Animation.DoubleAnimation
            $anim.From = 1.0
            $anim.To = 0.45
            $anim.Duration = [System.Windows.Duration]::new([TimeSpan]::FromMilliseconds(800))
            $anim.AutoReverse = $true
            $anim.RepeatBehavior = [System.Windows.Media.Animation.RepeatBehavior]::Forever
            $node.Dot.BeginAnimation([System.Windows.UIElement]::OpacityProperty, $anim)
            if ($PSBoundParameters.ContainsKey('Detail')) {
                $node.Detail.Text = $Detail
                $node.Detail.Visibility = if ($Detail) { 'Visible' } else { 'Collapsed' }
            }
            $node.Duration.Text = ''
        }
        'done' {
            # Stop the pulse animation
            $node.Dot.BeginAnimation([System.Windows.UIElement]::OpacityProperty, $null)
            $node.Dot.Opacity = 1.0
            $node.Dot.Background = New-Object System.Windows.Media.SolidColorBrush $green
            $node.Dot.BorderBrush = New-Object System.Windows.Media.SolidColorBrush $green
            $node.Dot.Effect = $null
            $node.Glyph.Text = [char]0x2713   # check
            $node.Glyph.Foreground = New-Object System.Windows.Media.SolidColorBrush $white
            $node.Label.Foreground = New-Object System.Windows.Media.SolidColorBrush $ink
            $node.Label.FontWeight = 'Normal'
            $node.Detail.Visibility = 'Collapsed'
            if ($Script:TaskTimers.ContainsKey($Key)) {
                $elapsed = (Get-Date) - $Script:TaskTimers[$Key]
                $node.Duration.Text = if ($elapsed.TotalMinutes -ge 1) {
                    "{0:N1}m" -f $elapsed.TotalMinutes
                } else {
                    "{0:N1}s" -f $elapsed.TotalSeconds
                }
            }
            $node.Duration.Foreground = New-Object System.Windows.Media.SolidColorBrush ([System.Windows.Media.ColorConverter]::ConvertFromString('#FFB0B8C8'))
        }
        'skip' {
            $node.Dot.BeginAnimation([System.Windows.UIElement]::OpacityProperty, $null)
            $node.Dot.Opacity = 1.0
            $node.Dot.Background = New-Object System.Windows.Media.SolidColorBrush $white
            $node.Dot.BorderBrush = New-Object System.Windows.Media.SolidColorBrush $muted
            $node.Dot.Effect = $null
            $node.Glyph.Text = [char]0x2212   # minus
            $node.Glyph.Foreground = New-Object System.Windows.Media.SolidColorBrush $muted
            $node.Label.Foreground = New-Object System.Windows.Media.SolidColorBrush $muted
            $node.Label.FontWeight = 'Normal'
            $node.Detail.Visibility = 'Collapsed'
            $node.Duration.Text = 'skipped'
            $node.Duration.Foreground = New-Object System.Windows.Media.SolidColorBrush $muted
        }
        'error' {
            $node.Dot.BeginAnimation([System.Windows.UIElement]::OpacityProperty, $null)
            $node.Dot.Opacity = 1.0
            $node.Dot.Background = New-Object System.Windows.Media.SolidColorBrush $red
            $node.Dot.BorderBrush = New-Object System.Windows.Media.SolidColorBrush $red
            $node.Dot.Effect = $null
            $node.Glyph.Text = [char]0x00D7   # multiplication sign / x
            $node.Glyph.Foreground = New-Object System.Windows.Media.SolidColorBrush $white
            $node.Label.Foreground = New-Object System.Windows.Media.SolidColorBrush $red
            $node.Label.FontWeight = 'SemiBold'
            if ($PSBoundParameters.ContainsKey('Detail') -and $Detail) {
                $node.Detail.Text = $Detail
                $node.Detail.Visibility = 'Visible'
                $node.Detail.Foreground = New-Object System.Windows.Media.SolidColorBrush $red
            }
            $node.Duration.Text = 'failed'
            $node.Duration.Foreground = New-Object System.Windows.Media.SolidColorBrush $red
        }
    }

    # Update counter (completed = done + skip)
    $completed = 0
    foreach ($k in $Script:TaskNodes.Keys) {
        $g = $Script:TaskNodes[$k].Glyph.Text
        if ($g -in @([char]0x2713, [char]0x2212, [char]0x00D7)) { $completed++ }
    }
    (Get-Element 'TaskCounter').Text = ("{0} / {1}" -f $completed, $Script:PipelineTasks.Count)

    [System.Windows.Forms.Application]::DoEvents()
}

# Adds a detail line under a running task (e.g. "downloading 1.2 / 4.0 MB").
# This is the new lightweight replacement for Add-SetupLog UI-side.
function Set-TaskDetail {
    param([string]$Key, [string]$Detail)
    $node = $Script:TaskNodes[$Key]
    if (-not $node) { return }
    $node.Detail.Text = $Detail
    $node.Detail.Visibility = if ($Detail) { 'Visible' } else { 'Collapsed' }
    [System.Windows.Forms.Application]::DoEvents()
}

# Add-SetupLog now only writes to the log file. The UI is driven by Update-Task /
# Set-TaskDetail. This preserves full debug history without the ugly terminal.
function Add-SetupLog {
    param([string]$Message)
    $stamp = (Get-Date -Format 'HH:mm:ss')
    Add-Content -Path $Script:SetupLogPath -Value "[$stamp] $Message"
}

function Set-SetupProgress {
    param([int]$Percent, [string]$Text)
    (Get-Element 'SetupProgress').Value = $Percent
    if ($Text) { (Get-Element 'SetupProgressText').Text = $Text.ToLower() }
    [System.Windows.Forms.Application]::DoEvents()
}

# Wire up events
(Get-Element 'NextBtn').Add_Click({
    Show-Page ($Script:CurrentPage + 1)
})

(Get-Element 'BackBtn').Add_Click({
    Show-Page ($Script:CurrentPage - 1)
})

(Get-Element 'CancelBtn').Add_Click({
    $window.Close()
})

(Get-Element 'MatlabBrowseBtn').Add_Click({
    $dlg = New-Object System.Windows.Forms.OpenFileDialog
    $dlg.Filter = 'matlab.exe|matlab.exe'
    $dlg.Title = 'Seleziona matlab.exe'
    if ($dlg.ShowDialog() -eq 'OK') {
        (Get-Element 'MatlabPathBox').Text = $dlg.FileName
    }
})

(Get-Element 'MatlabPathBox').Add_TextChanged({ Update-MatlabValidation })

(Get-Element 'OpenMathworksBtn').Add_Click({
    Start-Process 'https://www.mathworks.com/downloads/'
})

(Get-Element 'SnapBrowseBtn').Add_Click({
    $dlg = New-Object System.Windows.Forms.OpenFileDialog
    $dlg.Filter = 'gpt.exe|gpt.exe'
    $dlg.Title = 'Seleziona gpt.exe (in SNAP\bin\)'
    if ($dlg.ShowDialog() -eq 'OK') {
        (Get-Element 'SnapPathBox').Text = $dlg.FileName
    }
})

(Get-Element 'SnapPathBox').Add_TextChanged({ Update-SnapValidation })

(Get-Element 'InstallSnapBtn').Add_Click({
    (Get-Element 'InstallSnapBtn').IsEnabled = $false
    (Get-Element 'InstallSnapBtn').Content = 'Working...'
    # Visible feedback: a determinate bar during download, then indeterminate
    # marquee during install4j silent execution (install4j has no progress
    # API). The bar stays visible throughout so the user never sees a
    # frozen-looking "100%" state.
    (Get-Element 'SnapProgress').IsIndeterminate = $false
    (Get-Element 'SnapProgress').Value = 0
    (Get-Element 'SnapProgress').Visibility = 'Visible'
    (Get-Element 'SnapProgressText').Visibility = 'Visible'
    (Get-Element 'SnapProgressText').Text = 'preparing...'
    try {
        Invoke-SnapInstaller -InstallerPath $Script:BundledSnapPath `
            -StatusCallback {
                param($m)
                (Get-Element 'InstallSnapBtn').Content = $m
                (Get-Element 'SnapProgressText').Text = $m.ToLower()
            } `
            -ProgressCallback {
                param($pct)
                (Get-Element 'SnapProgress').Value = $pct
                (Get-Element 'SnapProgressText').Text = "downloading from esa $pct%"
            } `
            -InstallStartCallback {
                # Switch to indeterminate marquee — install4j silent does
                # not emit progress, but the user must see motion to know
                # the wizard is still working.
                (Get-Element 'SnapProgress').IsIndeterminate = $true
                (Get-Element 'SnapProgressText').Text = 'installing snap (UAC required, ~3-5 min)...'
                [System.Windows.Forms.Application]::DoEvents()
            } | Out-Null
        Start-Sleep -Seconds 2
        (Get-Element 'SnapProgress').IsIndeterminate = $false
        $found = Find-Snap
        if ($found) {
            (Get-Element 'SnapPathBox').Text = $found
            (Get-Element 'SnapStatus').Text = "[OK] SNAP installed successfully."
            (Get-Element 'SnapStatus').Foreground = '#FF008000'
            (Get-Element 'SnapInstallHint').Visibility = 'Collapsed'
        } else {
            (Get-Element 'InstallSnapBtn').Content = 'Retry detection'
            (Get-Element 'InstallSnapBtn').IsEnabled = $true
            (Get-Element 'SnapProgress').Visibility = 'Collapsed'
            (Get-Element 'SnapProgressText').Visibility = 'Collapsed'
            [System.Windows.MessageBox]::Show('SNAP installed but gpt.exe not detected. Provide the path manually using Browse.', 'SNAP detection', 'OK', 'Warning') | Out-Null
        }
    } catch {
        (Get-Element 'SnapProgress').IsIndeterminate = $false
        [System.Windows.MessageBox]::Show("Error during SNAP install:`n$($_.Exception.Message)", 'Error', 'OK', 'Error') | Out-Null
        (Get-Element 'InstallSnapBtn').Content = 'Install SNAP now'
        (Get-Element 'InstallSnapBtn').IsEnabled = $true
        (Get-Element 'SnapProgress').Visibility = 'Collapsed'
        (Get-Element 'SnapProgressText').Visibility = 'Collapsed'
    }
})

(Get-Element 'InstallPythonBtn').Add_Click({
    (Get-Element 'InstallPythonBtn').IsEnabled = $false
    (Get-Element 'PythonProgress').Visibility = 'Visible'
    (Get-Element 'PythonProgressText').Visibility = 'Visible'
    try {
        $statusCb = {
            param($msg)
            (Get-Element 'PythonProgressText').Text = $msg
            [System.Windows.Forms.Application]::DoEvents()
        }
        $progressCb = {
            param($pct)
            (Get-Element 'PythonProgress').Value = $pct
            [System.Windows.Forms.Application]::DoEvents()
        }
        $result = Install-PythonSilent -ProgressCallback $progressCb -StatusCallback $statusCb
        Install-PythonPackages -PythonExe $result.Exe -StatusCallback $statusCb
        $Script:State.PythonExe = $result.Exe
        $Script:State.PythonVersion = $result.Version
        (Get-Element 'PythonStatus').Text = "[OK] Python $($result.Version) installed."
        (Get-Element 'PythonStatus').Foreground = '#FF008000'
        (Get-Element 'PythonPath').Text = $result.Exe
        (Get-Element 'PythonInstallPanel').Visibility = 'Collapsed'
        (Get-Element 'NextBtn').IsEnabled = $true
    } catch {
        [System.Windows.MessageBox]::Show("Error during Python install:`n$($_.Exception.Message)", 'Error', 'OK', 'Error') | Out-Null
        (Get-Element 'InstallPythonBtn').IsEnabled = $true
    }
})

(Get-Element 'DestBrowseBtn').Add_Click({
    $dlg = New-Object System.Windows.Forms.FolderBrowserDialog
    $dlg.Description = 'Cartella di destinazione per PHASE'
    $dlg.SelectedPath = (Get-Element 'DestPathBox').Text
    if ($dlg.ShowDialog() -eq 'OK') {
        (Get-Element 'DestPathBox').Text = $dlg.SelectedPath
    }
})

(Get-Element 'DestPathBox').Add_TextChanged({ Update-DestValidation })

(Get-Element 'StartSetupBtn').Add_Click({
    (Get-Element 'StartSetupBtn').IsEnabled = $false
    (Get-Element 'CancelBtn').IsEnabled = $false
    try {
        Invoke-FullSetup
        Show-Page 7
    } catch {
        Add-SetupLog "ERRORE FATALE: $($_.Exception.Message)"
        [System.Windows.MessageBox]::Show("Installation failed:`n$($_.Exception.Message)`n`nLog: $Script:SetupLogPath", 'Error', 'OK', 'Error') | Out-Null
        (Get-Element 'CancelBtn').IsEnabled = $true
    }
})

(Get-Element 'OpenFolderBtn').Add_Click({
    Start-Process explorer.exe -ArgumentList (Join-Path $Script:State.InstallDir 'PHASE')
})

(Get-Element 'OpenLogBtn').Add_Click({
    Start-Process notepad.exe -ArgumentList $Script:SetupLogPath
})

# -----------------------------------------------------------------------------
# Create in the install root the .lnk shortcuts to the three .mlapp apps + a
# README, so the user launches the app from the main folder without browsing
# into engine\. The shortcuts point at the .mlapp files inside engine\PHASE\
# and set WorkingDirectory to the .mlapp folder (same cwd as a direct double
# click, so the relative .mat auto-load keeps working).
# -----------------------------------------------------------------------------
function New-PhaseLauncherShortcuts {
    param(
        [Parameter(Mandatory)] [string]$InstallDir,
        [Parameter(Mandatory)] [string]$PhaseDir,
        [scriptblock]$StatusCallback = { param($m) }
    )

    $apps = @(
        @{ Name = 'PHASE Preprocessing'; Target = (Join-Path $PhaseDir 'PHASE_Preprocessing.mlapp') }
        @{ Name = 'PHASE StaMPS';        Target = (Join-Path $PhaseDir 'PHASE_Preprocessing\PHASE_StaMPS.mlapp') }
        @{ Name = 'PHASE model';         Target = (Join-Path $PhaseDir 'PHASE_model.mlapp') }
    )

    $wsh = New-Object -ComObject WScript.Shell
    try {
        foreach ($a in $apps) {
            if (-not (Test-Path $a.Target)) {
                & $StatusCallback "[!] Shortcut skipped, target missing: $($a.Target)"
                continue
            }
            $lnkPath = Join-Path $InstallDir ($a.Name + '.lnk')
            $sc = $wsh.CreateShortcut($lnkPath)
            $sc.TargetPath = $a.Target
            $sc.WorkingDirectory = (Split-Path -Parent $a.Target)
            $sc.Description = "Open $($a.Name) in MATLAB App Designer"
            $sc.Save()
            & $StatusCallback "[OK] Shortcut: $lnkPath"
        }
    } finally {
        [System.Runtime.InteropServices.Marshal]::ReleaseComObject($wsh) | Out-Null
    }

    $readme = Join-Path $InstallDir 'README.txt'
    $readmeText = @"
PHASE - InSAR PSI suite
=======================

To START the application, double-click one of these shortcuts:

  - "PHASE Preprocessing.lnk"  ->  SNAP data preparation (module 1)
  - "PHASE StaMPS.lnk"         ->  StaMPS processing + export (module 2)
  - "PHASE model.lnk"          ->  modelling (module 3)

DATA INPUT
----------
  - Sentinel-1: retrieve images with the download script generated by
    the preprocessing module (Alaska SAR Facility).
  - Cosmo-SkyMed: in "PHASE Preprocessing" pick Cosmo-SkyMed and use the
    "Images" tab to import your .h5 files (they are copied into the
    slaves folder automatically).

GACOS atmospheric correction (optional)
---------------------------------------
  When you enable GACOS in "PHASE StaMPS", a window shows the request
  parameters (UTC, bounding box, dates) ready to paste into www.gacos.net.
  Select "Binary grid" as the file type, download the .tar.gz archives
  into the GACOS folder, then press Continue.

Do NOT move, rename or delete the "engine" folder: it holds PHASE's
code, native binaries and configuration. Moving it or deleting its
contents will prevent the application from starting.
"@
    Set-Content -Path $readme -Value $readmeText -Encoding UTF8
    & $StatusCallback "[OK] README: $readme"
}

# -----------------------------------------------------------------------------
# Main setup orchestration (runs in Page 6 when StartSetupBtn is clicked)
# -----------------------------------------------------------------------------

function Invoke-FullSetup {
    Add-SetupLog "=== PHASE Installer ==="
    Add-SetupLog "Log file: $Script:SetupLogPath"
    Add-SetupLog "MATLAB:  $($Script:State.MatlabExe)"
    Add-SetupLog "SNAP:    $($Script:State.SnapGpt)"
    Add-SetupLog "Python:  $($Script:State.PythonExe) ($($Script:State.PythonVersion))"
    Add-SetupLog "Dest:    $(Join-Path $Script:State.InstallDir 'PHASE')"

    # Task 1: git
    Set-SetupProgress 5 'detecting git'
    Update-Task -Key 'git' -Status 'running' -Detail 'looking for git on PATH...'
    $git = Find-Git
    if (-not $git) {
        Update-Task -Key 'git' -Status 'running' -Detail 'downloading Portable Git (~60 MB)...'
        Add-SetupLog "git not found on system. Downloading Portable Git..."
        try {
            $git = Install-PortableGit `
                -StatusCallback { param($m) Add-SetupLog $m; Set-TaskDetail -Key 'git' -Detail $m } `
                -ProgressCallback { param($pct) Set-SetupProgress $pct "downloading portable git $pct%"; Set-TaskDetail -Key 'git' -Detail "downloading $pct%" }
        } catch {
            Update-Task -Key 'git' -Status 'error' -Detail $_.Exception.Message
            throw "Unable to install Portable Git: $($_.Exception.Message). Install manually from git-scm.com/download/win and re-run the installer."
        }
    }
    $Script:State.GitExe = $git
    Update-Task -Key 'git' -Status 'done'
    Add-SetupLog "[OK] git: $git"

    # Always install into a dedicated "PHASE" folder inside the chosen
    # destination, so the picked folder isn't filled directly. This is the
    # user-facing app folder: it holds the .lnk shortcuts + README at the top.
    $appDir = Join-Path $Script:State.InstallDir 'PHASE'
    if (-not (Test-Path $appDir)) {
        New-Item -ItemType Directory -Path $appDir -Force | Out-Null
        Add-SetupLog "[OK] Created $appDir"
    }
    # Everything (repos + binaries + config) lives under PHASE\engine\ so the app
    # folder stays clean and only shows the .lnk shortcuts + README. The PHASE
    # repo is cloned directly into engine\ (no extra nesting), with StaMPS and
    # TRAIN as engine\StaMPS / engine\TRAIN. All downstream paths derive from
    # $phaseDir, so savepath / input_StaMPS.mat / project.conf follow along.
    $engineDir = Join-Path $appDir 'engine'
    $phaseDir = $engineDir                       # PHASE repo == engine\
    $stampsDir = Join-Path $phaseDir 'StaMPS'    # engine\StaMPS
    $trainDir = Join-Path $phaseDir 'TRAIN'      # engine\TRAIN

    # Task 2: clone PHASE
    Set-SetupProgress 15 'cloning phase'
    Update-Task -Key 'clone-phase' -Status 'running' -Detail $Script:PhaseRepo
    Invoke-GitClone -GitExe $git -Repo $Script:PhaseRepo -Branch $Script:PhaseBranch `
        -Destination $phaseDir -StatusCallback { param($m) Add-SetupLog $m; Set-TaskDetail -Key 'clone-phase' -Detail $m }
    Update-Task -Key 'clone-phase' -Status 'done'
    Add-SetupLog "[OK] PHASE cloned at $phaseDir"

    # Task 3: clone StaMPS
    Set-SetupProgress 30 'cloning stamps'
    Update-Task -Key 'clone-stamps' -Status 'running' -Detail $Script:StampsRepo
    Invoke-GitClone -GitExe $git -Repo $Script:StampsRepo -Branch $Script:StampsBranch `
        -Destination $stampsDir -StatusCallback { param($m) Add-SetupLog $m; Set-TaskDetail -Key 'clone-stamps' -Detail $m }
    Update-Task -Key 'clone-stamps' -Status 'done'
    Add-SetupLog "[OK] StaMPS cloned at $stampsDir"

    # Task 4: clone TRAIN
    Set-SetupProgress 45 'cloning train'
    Update-Task -Key 'clone-train' -Status 'running' -Detail $Script:TrainRepo
    Invoke-GitClone -GitExe $git -Repo $Script:TrainRepo -Branch $Script:TrainBranch `
        -Destination $trainDir -StatusCallback { param($m) Add-SetupLog $m; Set-TaskDetail -Key 'clone-train' -Detail $m }
    Update-Task -Key 'clone-train' -Status 'done'
    Add-SetupLog "[OK] TRAIN cloned at $trainDir"

    # Task 5: StaMPS native binaries (sempre da release precompilata - non
    # tentiamo piu' la build da sorgente con install-windows.ps1)
    Set-SetupProgress 60 'downloading native binaries'
    Update-Task -Key 'stamps-bin' -Status 'running' -Detail 'fetching stamps-win64-binaries.zip (~4 MB)...'
    $binOk = Invoke-StampsBinariesDownload -StampsRoot $stampsDir -StatusCallback { param($m) Add-SetupLog $m; Set-TaskDetail -Key 'stamps-bin' -Detail $m }
    if ($binOk) {
        Update-Task -Key 'stamps-bin' -Status 'done'
        Add-SetupLog "[OK] 7 StaMPS binaries ready in $stampsDir\bin"
    } else {
        Update-Task -Key 'stamps-bin' -Status 'error' -Detail 'download failed - mt_prep_snap will not work'
        Add-SetupLog "[!] StaMPS binaries download failed - mt_prep_snap will not work."
    }

    # Task 6: GMT
    Set-SetupProgress 75 'installing gmt'
    Update-Task -Key 'gmt' -Status 'running' -Detail 'checking GMT installation...'
    try {
        $gmtPath = Install-GmtSilent `
            -StatusCallback { param($m) Add-SetupLog $m; Set-TaskDetail -Key 'gmt' -Detail $m } `
            -ProgressCallback { param($pct) Set-SetupProgress $pct "downloading gmt $pct%"; Set-TaskDetail -Key 'gmt' -Detail "downloading $pct%" }
        Update-Task -Key 'gmt' -Status 'done'
        Add-SetupLog "[OK] GMT ready: $gmtPath"
    } catch {
        Update-Task -Key 'gmt' -Status 'skip' -Detail 'optional - only needed for a_gacos'
        Add-SetupLog "[!] GMT install failed: $($_.Exception.Message). Only required for tropo_method=a_gacos."
    }

    # Task 7: env configuration
    Set-SetupProgress 82 'configuring environment'
    Update-Task -Key 'env' -Status 'running' -Detail 'writing env vars and config files...'
    Set-MatlabEnvVar -MatlabExe $Script:State.MatlabExe
    Add-SetupLog "[OK] MATLAB_EXE set (user env var)"
    Set-PhasePythonConfig -PythonExe $Script:State.PythonExe
    Add-SetupLog "[OK] %APPDATA%\PHASE\python.txt written"
    Write-ProjectConfTemplate -InstallDir $phaseDir -SnapGpt $Script:State.SnapGpt
    Add-SetupLog "[OK] project.conf.template written at $phaseDir"
    Update-Task -Key 'env' -Status 'done'

    # Task 8: MATLAB savepath + .mat files
    Set-SetupProgress 90 'matlab savepath + .mat files'
    Update-Task -Key 'matlab' -Status 'running' -Detail 'launching matlab -batch (30-60s startup)...'
    $savepathResult = Invoke-MatlabSavePath -MatlabExe $Script:State.MatlabExe `
        -StampsRoot $stampsDir -TrainRoot $trainDir -PhaseRoot $phaseDir `
        -PythonExe $Script:State.PythonExe -SnapGpt $Script:State.SnapGpt `
        -StatusCallback { param($m) Add-SetupLog $m; Set-TaskDetail -Key 'matlab' -Detail $m }
    if ($savepathResult.Success) {
        Update-Task -Key 'matlab' -Status 'done'
        Add-SetupLog "[OK] MATLAB savepath OK (StaMPS + TRAIN added permanently)"
        Add-SetupLog "[OK] input_StaMPS.mat pre-populated (installation_folder + project_path)"
        Add-SetupLog "[OK] input_preprocessing.mat pre-populated (python + gptbin_path)"
    } else {
        Update-Task -Key 'matlab' -Status 'skip' -Detail "exit $($savepathResult.ExitCode) - run addpath/savepath manually"
        Add-SetupLog "[!] MATLAB savepath not confirmed (exit code $($savepathResult.ExitCode))."
        if ($savepathResult.Output) {
            Add-SetupLog "MATLAB output:"
            foreach ($line in ($savepathResult.Output -split "`n")) {
                if ($line.Trim()) { Add-SetupLog "    $line" }
            }
        }
        Add-SetupLog ""
        Add-SetupLog "Open MATLAB and run these commands manually:"
        Add-SetupLog "    addpath(genpath('$($stampsDir.Replace('\','/'))/matlab')); savepath"
        if ($trainDir) {
            Add-SetupLog "    addpath(genpath('$($trainDir.Replace('\','/'))/matlab')); savepath"
        }
    }

    # Task 9: patch .mlapp files
    Set-SetupProgress 95 'patching mlapp files'
    Update-Task -Key 'patch' -Status 'running' -Detail 'injecting startupFcn auto-load patches...'
    Get-Process matlab -ErrorAction SilentlyContinue | ForEach-Object {
        try { $_ | Stop-Process -Force; Add-SetupLog "MATLAB closed (PID $($_.Id)) to avoid stale class cache" } catch {}
    }
    Start-Sleep -Seconds 2
    Get-Process matlab -ErrorAction SilentlyContinue | ForEach-Object {
        try { $_ | Stop-Process -Force; Add-SetupLog "MATLAB chiuso (PID $($_.Id)) per evitare cache stale" } catch {}
    }
    Start-Sleep -Seconds 2
    # Le patch .mlapp sono SOLO comodita' (auto-load dei default all'apertura
    # dell'app). A questo punto tutta l'installazione vera (repo, binari, env
    # var, savepath, .mat) e' gia' completata. Un fallimento qui - p.es. un
    # anchor che non combacia perche' upstream ha cambiato lo startupFcn - NON
    # deve abortire l'intera installazione: degradiamo a warning e proseguiamo.
    try {
    $stampsMlapp = Join-Path $phaseDir 'PHASE_Preprocessing\PHASE_StaMPS.mlapp'
    [void](Invoke-MlappAutoLoadPatch -MlappPath $stampsMlapp `
        -MatFileRelative './input_StaMPS.mat' `
        -StatusCallback { param($m) Add-SetupLog $m })

    # PHASE_Preprocessing.mlapp: fixa 2 problemi all'avvio dell'app:
    #   1. Entrambi i pannelli (Sentinel1 + CosmoSkyMed) sono invisibili by
    #      default - l'utente deve muovere lo switch per vedere qualcosa.
    #      Fix: attiviamo Sentinel1 (default sensato).
    #   2. Solo i 2 path di INSTALLAZIONE (Python + GPT) vanno precompilati;
    #      gli altri 21 campi sono dataset-specifici. Anziche' chiamare il
    #      LoadButton (che fallisce silenziosamente su un campo intermedio),
    #      settiamo direttamente python_SEN/gptbin_path_SEN dal .mat se esiste.
    $prepInject = @"


            % AUTO-CONFIG (PHASE installer): default pannello SEN + path Python/GPT
            try
                % 1. Default Sentinel1 panel visibile
                if isprop(app, 'ConstellationSwitch') && isvalid(app.ConstellationSwitch)
                    app.ConstellationSwitch.Value = 'Sentinel1';
                end
                if isprop(app, 'Sentinel1Panel') && isvalid(app.Sentinel1Panel)
                    app.Sentinel1Panel.Visible = 'on';
                end
                if isprop(app, 'CosmoSkyMedPanel') && isvalid(app.CosmoSkyMedPanel)
                    app.CosmoSkyMedPanel.Visible = 'off';
                end
                app.constellation = 'SEN';
                drawnow;

                % 2. Set diretto dei 2 path da input_preprocessing.mat (Python+GPT)
                if exist('./PHASE_Preprocessing/input_preprocessing.mat', 'file') == 2
                    cfg = load('./PHASE_Preprocessing/input_preprocessing.mat');
                    if isfield(cfg, 'python')
                        % SEN side
                        if isprop(app, 'CustomPythonEnvironmentEditField') && isvalid(app.CustomPythonEnvironmentEditField)
                            app.CustomPythonEnvironmentEditField.Value = cfg.python;
                            app.CustomPythonEnvironmentEditField.Visible = 'on';
                        end
                        if isprop(app, 'PythonEnvironmentDropDown') && isvalid(app.PythonEnvironmentDropDown)
                            app.PythonEnvironmentDropDown.Value = 'Other';
                        end
                        if isprop(app, 'PythonEnvironmentLabel') && isvalid(app.PythonEnvironmentLabel)
                            app.PythonEnvironmentLabel.Visible = 'on';
                        end
                        app.python_SEN = cfg.python;
                        % CSK side (stesso python_exe)
                        if isprop(app, 'CustomPythonEnvironmentEditField_2') && isvalid(app.CustomPythonEnvironmentEditField_2)
                            app.CustomPythonEnvironmentEditField_2.Value = cfg.python;
                            app.CustomPythonEnvironmentEditField_2.Visible = 'on';
                        end
                        if isprop(app, 'PythonEnvironmentDropDown_2') && isvalid(app.PythonEnvironmentDropDown_2)
                            app.PythonEnvironmentDropDown_2.Value = 'Other';
                        end
                        app.python_CSK = cfg.python;
                    end
                    if isfield(cfg, 'gptbin_path')
                        if isprop(app, 'PathEditField') && isvalid(app.PathEditField)
                            app.PathEditField.Value = cfg.gptbin_path;
                        end
                        app.gptbin_path_SEN = cfg.gptbin_path;
                        if isprop(app, 'PathEditField_2') && isvalid(app.PathEditField_2)
                            app.PathEditField_2.Value = cfg.gptbin_path;
                        end
                        app.gptbin_path_CSK = cfg.gptbin_path;
                    end
                end
            catch
                % non blocca lo startup se qualcosa fallisce
            end
"@
    $prepMlapp = Join-Path $phaseDir 'PHASE_Preprocessing.mlapp'
    # Anchor specifico: l'ultima delle 4 righe di "Initially hide" (linea 479
    # del document.xml originale). Inserire qui assicura che la nostra
    # visibility='on' sui field Python non venga sovrascritta dalle 'off'
    # subito sotto la riga di anchor di default.
    [void](Invoke-MlappAutoLoadPatch -MlappPath $prepMlapp `
        -Anchor "app.CustomPythonEnvironmentEditField_2.Visible = 'off';" `
        -InjectBlock $prepInject `
        -StatusCallback { param($m) Add-SetupLog $m })

    # PHASE_model.mlapp ha gia' un auto-load nel suo startupFcn (legge config
    # da input_model.mat se esiste). Ma noi NON generiamo input_model.mat
    # perche' il config struct ha 60+ campi. Quindi patchiamo lo startupFcn
    # per settare SOLO app.pythonPath al python configurato, lasciando intatto
    # il resto del flusso (compreso il load esistente se il .mat sara' creato
    # poi dall'utente con Save).
    $pyForModel = $Script:State.PythonExe.Replace('\','/')
    $modelInject = @"


            % AUTO-CONFIG (PHASE installer): default pythonPath se input_model.mat
            % non esiste ancora. Sovrascritto dal config.pythonPath del load
            % successivo (linea ~341 dello startupFcn).
            try
                appDir_phaseinstaller = fileparts(mfilename('fullpath'));
                if exist(fullfile(appDir_phaseinstaller, 'input_model.mat'), 'file') ~= 2
                    app.pythonPath = '$pyForModel';
                    if isprop(app, 'pythoninstallationpathEditField') && isvalid(app.pythoninstallationpathEditField)
                        app.pythoninstallationpathEditField.Value = '$pyForModel';
                    end
                end
            catch
            end
"@
    $modelMlapp = Join-Path $phaseDir 'PHASE_model.mlapp'
    # Anchor: prima riga dello startupFcn (IndexOf prende la prima delle 2
    # occorrenze, che e' quella in startupFcn - non quella nel callback Run).
    # NB: deve combaciare LETTERALMENTE col codice MATLAB embedded nel
    # document.xml. Upstream usa fullfile(pwd,...) - non './MatlabFunctions/'.
    [void](Invoke-MlappAutoLoadPatch -MlappPath $modelMlapp `
        -Anchor "addpath(fullfile(pwd, 'MatlabFunctions'));" `
        -InjectBlock $modelInject `
        -StatusCallback { param($m) Add-SetupLog $m })

        Update-Task -Key 'patch' -Status 'done'
    } catch {
        # Patch non riuscita (tipicamente anchor non trovato per drift upstream
        # dei .mlapp). L'installazione resta valida: gli app si aprono comunque,
        # l'utente caricara' i default manualmente dal tab Save/Load.
        Add-SetupLog "WARNING: patch .mlapp non applicata ($($_.Exception.Message)) - installazione comunque valida, caricare i default manualmente."
        Update-Task -Key 'patch' -Status 'skip' -Detail 'anchor non trovato - apri il .mlapp e usa Load manualmente'
    }

    # Collegamenti + README nella cartella principale (engine\ resta nascosto
    # all'uso quotidiano). Non bloccante: se fallisce, i .mlapp restano comunque
    # apribili da engine\PHASE\.
    try {
        New-PhaseLauncherShortcuts -InstallDir $appDir -PhaseDir $phaseDir `
            -StatusCallback { param($m) Add-SetupLog $m }
    } catch {
        Add-SetupLog "[!] Could not create root shortcuts/README: $($_.Exception.Message)"
    }

    Set-SetupProgress 100 'all done'
    Add-SetupLog ""
    Add-SetupLog "=== Installation complete ==="
    Add-SetupLog "Launch the app from the shortcuts in $($appDir):"
    Add-SetupLog "  PHASE Preprocessing.lnk"
    Add-SetupLog "  PHASE StaMPS.lnk"
    Add-SetupLog "  PHASE model.lnk"
    Add-SetupLog "(the actual files live in $phaseDir - no need to open them by hand)"
}

# -----------------------------------------------------------------------------
# Show wizard
# -----------------------------------------------------------------------------

if ($DryRun) {
    Write-Host "Dry run: XAML parsed OK, $($pages.Count) pages registered."
    Write-Host "Detection probes:"
    Write-Host "  MATLAB: $(Find-Matlab)"
    Write-Host "  SNAP:   $(Find-Snap)"
    Write-Host "  Python: $((Find-Python | ConvertTo-Json -Compress))"
    Write-Host "  git:    $(Find-Git)"
    exit 0
}

Show-Page 1
[void]$window.ShowDialog()
