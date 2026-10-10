# Installer PHASE

> **Stato:** queste istruzioni per l'hub PHASE 7 descrivono il branch di sviluppo.
> La release pubblica 6.1.4 usa ancora i tre launcher precedenti. Non è ancora
> stata pubblicata una release PHASE 7 né validata una pipeline PSI completa su
> ciascun sistema operativo.

## Installazione dell'hub PHASE 7

Una sola installazione del codice serve più progetti. L'hub crea/apre cartelle
di progetto separate, con dati e output nel progetto; il motore MATLAB resta
nell'installazione. MATLAB e ESA SNAP devono essere già installati e licenziati
se necessario. Per il calcolo PSI servono inoltre runtime StaMPS e, se usato,
TRAIN funzionanti sul sistema: l'installer macOS/Linux può incorporare copie
**già preparate**, ma non le compila o ne certifica il funzionamento.

### macOS Apple Silicon e Linux

Richiede Python 3.10+ e `git` per il download da GitHub. Da un checkout di
sviluppo si può installare esattamente il codice locale (incluse modifiche non
ancora pubblicate) usando `--source`:

```bash
# macOS Apple Silicon
./installer/install-phase-macos.command --source "$PWD" \
  --matlab /Applications/MATLAB_R2026a.app/bin/matlab \
  --gpt /Applications/esa-snap/bin/gpt

# Linux: adattare i due percorsi alla propria installazione
./installer/install-phase-linux.sh --source "$PWD" \
  --matlab /path/to/matlab/bin/matlab --gpt /path/to/esa-snap/bin/gpt
```

Senza `--source`, lo script clona il branch `main` da GitHub, che deve già
contenere l'hub. Per vedere destinazioni e dipendenze senza scrivere file,
aggiungere `--dry-run`. `--stamps /path/to/StaMPS` e `--train /path/to/TRAIN`
incorporano runtime preparati; `--python` sceglie l'interprete con cui creare
l'ambiente isolato. `--skip-python-deps` è solo per sviluppo/test e non prepara
le librerie Python necessarie al processing.

La destinazione predefinita è `~/Library/Application Support/PHASE` su macOS e
`~/.local/share/PHASE` su Linux; `--prefix` permette un'altra cartella utente.
Lo script copia un runtime ripulito dai file di sviluppo, crea un ambiente
Python e installa `openpyxl`, `requests`, `asf_search`, `shapely`. Su macOS crea
`~/Applications/PHASE.app`; su Linux un launcher nel menu applicazioni e, se
libero, `~/.local/bin/phase`. In entrambi i casi l'avvio reale è
`<prefix>/launch-phase.sh`. In caso di aggiornamento, i componenti precedenti
sono conservati in `<prefix>/backups/`; le cartelle dei progetti non vengono
toccate. Se esiste già un collegamento `PHASE.app` o `phase.desktop` non gestito,
l'installer lo lascia intatto e indica il launcher diretto.

L'installer macOS rifiuta Apple Intel. Linux non è ancora stato testato con
un'installazione completa su una macchina Linux; il wrapper e la preparazione
dei file non equivalgono a una verifica end-to-end del processing.

#### Stato del porting PSI su macOS Apple Silicon

Il fork `pyccino/StaMPS` usato da PHASE è stato compilato su Apple Silicon:
sette binari arm64 e sette test CTest passano. SNAPHU, Triangle e GNU awk
sono stati compilati per arm64 e inclusi nel runtime di prova. Il launcher Unix imposta ora
`STAMPS`, `APS_toolbox`, MATLAB e i percorsi dei binari; l'installer sostituisce
le configurazioni upstream con percorsi locali reali. I comandi di preparazione
StaMPS falliscono esplicitamente se un prerequisito o il processing fallisce,
anche quando i progetti contengono spazi nel percorso.

Per preparare un runtime di prova prima di integrare tutto nel DMG:

```bash
brew install cmake
git clone https://github.com/pyccino/StaMPS.git /path/to/StaMPS-source
git clone https://github.com/pyccino/TRAIN.git /path/to/TRAIN-source
python3 installer/prepare-macos-runtime.py \
  --stamps-source /path/to/StaMPS-source \
  --train-source /path/to/TRAIN-source \
  --snaphu /path/to/arm64/snaphu \
  --triangle /path/to/arm64/triangle \
  --gawk /path/to/arm64/gawk \
  --output /path/to/phase-macos-runtime
./installer/install-phase-macos.command --source "$PWD" \
  --stamps /path/to/phase-macos-runtime/StaMPS \
  --train /path/to/phase-macos-runtime/TRAIN
```

SNAPHU, GNU awk e il programma Triangle di Shewchuk vanno ottenuti separatamente dalle
[pagine Stanford](https://web.stanford.edu/group/radar/softwareandlinks/sw/snaphu/)
e [CMU](https://www.cs.cmu.edu/~quake/triangle.html), e dal
[progetto GNU awk](https://www.gnu.org/software/gawk/), e compilati per arm64.
Non usare il pacchetto Homebrew `triangle`: è un programma omonimo diverso.
Questi componenti hanno condizioni di licenza da verificare prima di includerli in un
DMG pubblico. Per una **preview locale** che includa il runtime già preparato:

```bash
python3 installer/build-macos-dmg.py \
  --runtime /path/to/phase-macos-runtime \
  --output /path/to/PHASE-7-macos-arm64-preview.dmg
```

Senza `--runtime`, il DMG installa il solo motore PHASE e lo dichiara nel
wizard. Con `--runtime`, la GUI installa StaMPS/TRAIN insieme all'app, senza
chiedere all'utente di scegliere manualmente le loro cartelle. La preview non
è una release pubblica: mancano una prova della pipeline PSI su dati reali
macOS, una revisione delle licenze/distribuzione e firma/notarizzazione.

### Pacchetti grafici da distribuire

- **Windows:** compilare `install-phase.exe` su Windows con
  `compile-to-exe.ps1`. Il wizard scarica le dipendenze Windows che gestisce
  già e crea `PHASE.lnk` sia sul desktop sia nella cartella installata,
  usando `Logo_square.png` convertito in `PHASE.ico`. Il `.exe` non incorpora
  MATLAB o la licenza. Non è ancora stato compilato o provato su Windows per
  PHASE 7.
  Per una prova del branch prima della release, il workflow
  `Build PHASE 7 Windows test installer` produce un artefatto EXE con
  `codex/phase-stamps-beta` incorporato. Lo ZIP degli artefatti GitHub contiene
  un solo EXE; non serve passare argomenti a riga di comando all'utente.
- **macOS Apple Silicon:** costruire `PHASE-7-macos-arm64.dmg` con
  `python3 installer/build-macos-dmg.py --output /path/PHASE-7-macos-arm64.dmg`.
  Il DMG contiene un wizard nativo SwiftUI a sette passi, con palette e struttura
  visiva dell'installer Windows, e una copia del motore;
  crea poi `~/Applications/PHASE.app` con icona PHASE. MATLAB, SNAP e Python 3.10+
  con `venv` sono prerequisiti esterni. Per una distribuzione senza avvisi di
  Gatekeeper, usare `--sign-identity` e `--notary-profile` con credenziali
  Apple Developer ID configurate. Per compilare il wizard serve la toolchain
  Swift di Xcode; l'utente finale non deve installare Xcode. Il DMG costruito
  senza firma/notarizzazione è
  solo una preview locale.
- **Linux:** su Linux x86_64 o aarch64, con `appimagetool`, costruire
  `PHASE-7-linux.AppImage` tramite
  `python3 installer/build-linux-appimage.py --output /path/PHASE-7-linux.AppImage`.
  All'avvio mostra un'interfaccia con `zenity` o `kdialog` e installa il motore
  nel profilo utente. Richiede Python 3 con `venv`, MATLAB e SNAP; StaMPS/TRAIN
  vanno forniti come runtime già preparati. Builder e AppImage richiedono ancora
  una prova su Linux.

I pacchetti macOS/Linux includono il codice PHASE al momento della build, per
evitare che un aggiornamento successivo di `main` cambi ciò che installano.
Nessun pacchetto PHASE 7 è ancora allegato a una release pubblica.

Su macOS PHASE limita a ogni avvio di SNAP GPT la cache e il parallelismo in
base alla RAM e all'heap Java configurato, senza cambiare i parametri
scientifici: su un Mac da 8 GB con heap GPT da 5 GB, una configurazione
`-c 26G -q 8` viene eseguita come `-c 512M -q 2 -x` e il valore effettivo
compare nel Run monitor. Dopo un errore di memoria nello Step 3, mantenere i
prodotti degli Step 1–2 e ripartire da **First preprocessing step = 3**.
Il limite riduce la pressione sulla memoria, ma il completamento di una
pipeline su 8 GB dipende anche dalle dimensioni delle acquisizioni e dell'AOI.

### Aggiornamenti dall'hub

Il pulsante **Cerca update** interroga l'ultima release stabile GitHub di
`robimonti/PHASE`. Richiede un asset `phase7-engine.zip` della serie v7 e una
impronta SHA-256 negli asset della release. L'hub scarica, verifica e prepara
l'archivio; il launcher applica l'aggiornamento al successivo avvio, prima di
caricare MATLAB. La versione precedente viene conservata in `backups`. StaMPS,
TRAIN e i file di configurazione creati dall'installer restano nell'installazione.
I progetti esterni non vengono modificati. Per applicare l'update occorre
chiudere MATLAB e avviare PHASE dal collegamento installato, non digitare
`PHASE_Hub` in una sessione MATLAB già aperta.

Finché la prima release v7 con l'asset dedicato non è pubblicata, il pulsante
non proporrà aggiornamenti. Per preparare quell'asset da un checkout di release:

```bash
python3 installer/build-update-package.py --tag v7.0.0 \
  --output /path/outside/repository/phase7-engine.zip
```

Il builder richiede un checkout pulito esattamente sul tag indicato; `--dev-build`
serve solo per prove locali. L'archivio va allegato alla release GitHub con il medesimo tag. Il pacchetto
contiene solo il motore PHASE, non StaMPS, TRAIN, MATLAB, SNAP o dati di progetto.
L'aggiornamento non installa nuove dipendenze esterne: se una release ne richiede,
le istruzioni di release devono indicarle e l'installer va aggiornato.

### Windows

Il wizard PowerShell continua a rilevare/installare le dipendenze Windows già
supportate, ma ora propone `%LOCALAPPDATA%\Programs\PHASE` e crea un solo
collegamento `PHASE.lnk` per l'hub. Chi aggiorna un'installazione precedente
può ancora vedere i vecchi collegamenti: non vengono cancellati automaticamente.
Per provarlo dal sorgente, vedere la sezione seguente. La compilazione dell'EXE
e il test effettivo del wizard richiedono Windows; non sono stati eseguiti su
questo host macOS.

## Installer Windows storico e packaging

Wizard end-to-end (GUI WPF) del branch di sviluppo che installa PHASE 7 e le sue dipendenze su
Windows: MATLAB detection, SNAP install, Python 3.11+ silent install, clone di
PHASE/StaMPS/TRAIN, download verificato dei binari Triangle/snaphu, configurazione `MATLAB_EXE` +
`python.txt` + `savepath`. Per default clona il branch `main`; il branch può essere sovrascritto con
`-PhaseBranch`.

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

Produce `install-phase.exe` accanto allo script.

### 2. Bundle l'installer SNAP

Lo script cerca l'installer SNAP in `.\installers\esa-snap_sentinel_windows-13.0.0.exe`
accanto a sé. Per distribuirlo come pacchetto self-contained:

```powershell
# Layout finale del pacchetto:
phase-installer-v6.1.4\
├── install-phase.exe
└── installers\
    └── esa-snap_sentinel_windows-13.0.0.exe               # ~500 MB

# Comprimi:
Compress-Archive -Path phase-installer-v6.1.4 -DestinationPath phase-installer-v6.1.4.zip
```

L'utente finale estrae lo zip e fa doppio click su
`install-phase.exe`.

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
│   ├─ Invoke-StampsBinariesDownload scarica i 9 eseguibili Windows obbligatori
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
5. **Cartella destinazione** — default `%LOCALAPPDATA%\Programs\PHASE`.
   Validazione: scrivibile, no OneDrive (warning, non blocco), no caratteri
   non-ASCII.
6. **Installazione** — clona il branch `main` di PHASE + StaMPS + TRAIN, scarica e verifica i nove
   eseguibili StaMPS Windows (incluso `snaphu.exe`; un fallimento interrompe
   l'installazione), scrive `MATLAB_EXE` env var, scrive
   `%APPDATA%\PHASE\python.txt`, scrive `project.conf.template`, lancia
   `matlab.exe -batch` per addpath+savepath, rimuove dal runtime `legacy` e i
   file di sviluppo. Log live in console scrollabile.
7. **Fine** — riepilogo + bottoni "Apri cartella PHASE" e "Apri log".

La cartella `PHASE` contiene `PHASE.lnk` e il motore nella sottocartella
`engine`. I tre moduli sono sezioni dell'hub; i launcher standalone restano nel
motore per compatibilità durante la migrazione.

## Path configurati automaticamente

Dopo che l'installer ha finito, l'utente trova:

| Cosa | Dove | Valore |
|---|---|---|
| `MATLAB_EXE` env var (user scope) | `setx` registry | Path a `matlab.exe` |
| Python override per StaMPS | `%APPDATA%\PHASE\python.txt` | Path a `python.exe` (letto da `mt_prep_snap.bat:27`) |
| MATLAB path permanente (`pathdef.m`) | `matlab.exe -batch savepath` | `StaMPS\matlab` + `matlab_compat` + `TRAIN\matlab` |
| Template config dataset | `<dest>\PHASE\engine\project.conf.template` | `GPTBIN_PATH` precompilato + AOI placeholder |

L'utente avvia PHASE dal collegamento dell'hub. Il collegamento apre MATLAB,
aggiunge il motore al path ed esegue immediatamente la funzione standalone;
non apre il file nell'Editor.

## Caveat noti

1. **MATLAB non installabile automaticamente**: proprietario + licensing.
   L'installer fa solo detection + config.
2. **SNAP semi-interattivo**: l'installer ESA non ha modalità completamente
   silent senza response file pre-generato. Lo lanciamo standard, l'utente
   clicca Avanti×3 (~5 minuti).
3. **git**: se assente, il wizard installa Portable Git nel profilo utente.
4. **SmartScreen**: vedi sezione "Firma digitale" sopra.
5. **`matlab.exe -batch savepath`**: richiede licenza MATLAB già attivata.
   Se la licenza non è ancora stata accettata, il savepath fallisce con
   warning ma l'install procede; l'utente fa addpath/savepath manualmente
   alla prima apertura di MATLAB.
6. **Build dell'EXE**: PS2EXE richiede Windows PowerShell. Il sorgente
   `install-phase.ps1` può essere verificato nel repository, ma l'EXE finale
   va compilato su Windows con `compile-to-exe.ps1`.

## Disinstallazione

L'installer non scrive un uninstaller. Per pulire:

```powershell
# 1. Cancella cartella PHASE
Remove-Item -Recurse -Force "$env:LOCALAPPDATA\Programs\PHASE"

# 2. Rimuovi env var
[Environment]::SetEnvironmentVariable('MATLAB_EXE', $null, 'User')

# 3. Rimuovi config PHASE
Remove-Item -Recurse -Force "$env:APPDATA\PHASE"

# 4. (Opzionale) disinstalla Python 3.11.9 e SNAP da Pannello di Controllo
```
