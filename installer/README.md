# PHASE Windows Installer

Wizard end-to-end (GUI WPF) che installa PHASE e tutte le sue dipendenze su
Windows: MATLAB detection, SNAP install, Python 3.11+ silent install, clone di
PHASE/StaMPS/TRAIN, build di Triangle/snaphu, configurazione `MATLAB_EXE` +
`python.txt` + `savepath`.

## File

| File | Cosa è |
|---|---|
| `install-phase.ps1` | Sorgente PowerShell con WPF inline. Lanciabile direttamente. |
| `compile-to-exe.ps1` | Helper per compilare il `.ps1` in `.exe` via PS2EXE. |
| `README.md` | Questo file. |

## Uso (sorgente, sviluppo)

```powershell
# Una tantum, abilita gli script:
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned -Force

# Lancia il wizard:
powershell -ExecutionPolicy Bypass -File install-phase.ps1
```

Per debugging senza GUI:

```powershell
powershell -ExecutionPolicy Bypass -File install-phase.ps1 -DryRun
```

Stampa il risultato dei detector (MATLAB / SNAP / Python / git) ed esce con
codice 0.

## Uso (distribuzione)

### 1. Compila il `.ps1` in `.exe`

```powershell
# Una tantum:
Install-Module -Name ps2exe -Scope CurrentUser -Force

# Compila:
powershell -ExecutionPolicy Bypass -File compile-to-exe.ps1
```

Produce `install-phase.exe` (~3 MB) accanto allo script.

### 2. Bundle l'installer SNAP

Lo script cerca l'installer SNAP in `.\installers\esa-snap_sentinel_windows-13.0.0.exe`
accanto a sé. Per distribuirlo come pacchetto self-contained:

```powershell
# Layout finale del pacchetto:
phase-installer-v1.0.0\
├── install-phase.exe                                       # 3 MB
└── installers\
    └── esa-snap_sentinel_windows-13.0.0.exe               # ~500 MB

# Comprimi:
Compress-Archive -Path phase-installer-v1.0.0 -DestinationPath phase-installer-v1.0.0.zip
```

L'utente finale estrae lo zip e fa doppio click su `install-phase.exe`.

### 3. SmartScreen / firma digitale

L'`.exe` non firmato triggera **"Windows ha protetto il tuo PC"** al primo
lancio. L'utente clicca **Altre info → Esegui comunque** (una volta sola
per file, per utente).

Per evitarlo serve un certificato code-signing:
- A pagamento: ~€200/anno (Sectigo, DigiCert, Comodo).
- Gratuito per OSS: programma SignPath, lo stesso che StaMPS sta già usando
  (vedi `StaMPS/docs/SIGNPATH_STATUS.md`). Iter ~2-4 settimane.

## Architettura

```
install-phase.ps1
├─ Constants                 URL repo, branch, Python version, ecc.
├─ Detection helpers
│   ├─ Find-Matlab          Program Files glob + registry HKLM Mathworks + PATH
│   ├─ Find-Snap            Glob in C:\Program Files\esa-snap*\bin\gpt.exe + PATH
│   ├─ Find-Python          py launcher per minor 11..20, esclude \WindowsApps\
│   └─ Find-Git             Get-Command + path standard
├─ Action helpers
│   ├─ Get-RemoteFile       HTTP download con progress callback
│   ├─ Install-PythonSilent Lancia python-3.11.9-amd64.exe /quiet
│   ├─ Invoke-SnapInstaller Lancia esa-snap_sentinel_windows-13.0.0.exe (semi-interattivo)
│   ├─ Invoke-GitClone      git clone con branch e callback
│   ├─ Invoke-MatlabSavePath matlab.exe -batch addpath+savepath
│   ├─ Invoke-StampsInstall  Lancia StaMPS\install-windows.ps1
│   ├─ Set-MatlabEnvVar     setx MATLAB_EXE user scope
│   ├─ Set-PhasePythonConfig %APPDATA%\PHASE\python.txt
│   └─ Write-ProjectConfTemplate project.conf.template con GPTBIN_PATH
├─ WPF XAML                 7 pagine: Welcome, MATLAB, SNAP, Python, Dest, Setup, Finish
├─ Event handlers           Wire up dei click + validazione campi
└─ Invoke-FullSetup         Orchestratore (chiamato da Page 6 "Avvia installazione")
```

## Wizard step-by-step

1. **Welcome** — logo + intro.
2. **MATLAB** — auto-detect via Program Files glob, registry HKLM Mathworks,
   PATH. TextBox + "Sfoglia…". Avanti disabilitato finché path non valido.
   Se assente: link a mathworks.com (l'installer non può installare MATLAB
   perché proprietario).
3. **SNAP** — auto-detect via glob `C:\Program Files\esa-snap*\bin\gpt.exe`.
   Se assente e l'installer ESA è bundled: bottone "Installa SNAP ora" che
   lancia l'installer (semi-interattivo, l'utente clicca Avanti×3).
4. **Python** — auto-detect via `py -3.X` (X=11..20) escludendo `\WindowsApps\`.
   Se assente: download da python.org + silent install per-user (`/quiet
   InstallAllUsers=0 PrependPath=1`) con progress bar. Poi
   `pip install openpyxl requests asf_search shapely`.
5. **Cartella destinazione** — default `%USERPROFILE%\Desktop\PHASE`.
   Validazione: scrivibile, no OneDrive (warning, non blocco), no caratteri
   non-ASCII.
6. **Installazione** — clone PHASE + StaMPS + TRAIN, esegue
   `StaMPS\install-windows.ps1`, scrive `MATLAB_EXE` env var, scrive
   `%APPDATA%\PHASE\python.txt`, scrive `project.conf.template`, lancia
   `matlab.exe -batch` per addpath+savepath. Log live in console scrollabile.
7. **Fine** — riepilogo + bottoni "Apri cartella PHASE" e "Apri log".

## Path configurati automaticamente

Dopo che l'installer ha finito, l'utente trova:

| Cosa | Dove | Valore |
|---|---|---|
| `MATLAB_EXE` env var (user scope) | `setx` registry | Path a `matlab.exe` |
| Python override per StaMPS | `%APPDATA%\PHASE\python.txt` | Path a `python.exe` (letto da `mt_prep_snap.bat:27`) |
| MATLAB path permanente (`pathdef.m`) | `matlab.exe -batch savepath` | `StaMPS\matlab` + `matlab_compat` + `TRAIN\matlab` |
| Template config dataset | `<dest>\PHASE\project.conf.template` | `GPTBIN_PATH` precompilato + AOI placeholder |

L'utente apre uno qualsiasi dei `.mlapp` da MATLAB e tutto funziona. Per ogni
nuovo dataset deve solo copiare `project.conf.template` in `project.conf` e
riempire `MASTER` + bounding box AOI.

## Caveat noti

1. **MATLAB non installabile automaticamente**: proprietario + licensing.
   L'installer fa solo detection + config.
2. **SNAP semi-interattivo**: l'installer ESA non ha modalità completamente
   silent senza response file pre-generato. Lo lanciamo standard, l'utente
   clicca Avanti×3 (~5 minuti).
3. **git richiesto sul sistema**: per il clone. Se non c'è, l'installer
   stoppa con messaggio actionable. (TODO: download di Portable Git
   on-the-fly).
4. **SmartScreen**: vedi sezione "Firma digitale" sopra.
5. **`matlab.exe -batch savepath`**: richiede licenza MATLAB già attivata.
   Se la licenza non è ancora stata accettata, il savepath fallisce con
   warning ma l'install procede; l'utente fa addpath/savepath manualmente
   alla prima apertura di MATLAB.

## Disinstallazione

L'installer non scrive un uninstaller. Per pulire:

```powershell
# 1. Cancella cartella PHASE
Remove-Item -Recurse -Force "$env:USERPROFILE\Desktop\PHASE"

# 2. Rimuovi env var
[Environment]::SetEnvironmentVariable('MATLAB_EXE', $null, 'User')

# 3. Rimuovi config PHASE
Remove-Item -Recurse -Force "$env:APPDATA\PHASE"

# 4. (Opzionale) disinstalla Python 3.11.9 e SNAP da Pannello di Controllo
```
