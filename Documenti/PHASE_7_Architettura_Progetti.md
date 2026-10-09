# PHASE 7 — Contratto dei progetti (v1)

## Decisioni

- Una installazione PHASE contiene codice e risorse comuni. Ogni progetto contiene dati, configurazioni, elaborazioni, risultati e log. Il progetto non contiene una copia di PHASE.
- Piattaforme previste: Windows, Linux e macOS Apple Silicon. macOS Intel non fa parte della matrice di supporto.
- `phase-project.json` identifica un progetto. La versione dello schema è separata dalla versione dell'app.
- I file del progetto usano percorsi relativi alla sua cartella. I riferimenti a file esterni, per esempio quelli di un vecchio progetto, sono indicati come tali nel catalogo importazioni.
- La migrazione importa i prodotti scientifici finali. Non promette di riprendere un'elaborazione legacy interrotta.

## Layout v1

```text
Progetto/
  phase-project.json
  input/
    raw/
    aoi/
  processing/
    preprocessing/
    stamps/
  results/
    exports/
    model/
    figures/
    gis/
    reports/
  logs/
  imports/
    legacy-import.json       # presente dopo un'importazione
```

`phase-project.json` contiene `schemaVersion`, `id` (UUID), `name`, `createdAt` (UTC) e `layout`. La v1 non memorizza percorsi assoluti per le cartelle proprie del progetto: `phase_project.paths` li ricava dalla cartella che contiene il manifest. Il progetto può quindi essere spostato.

Il codice MATLAB usa `phase_project.installationRoot()` per le risorse dell'app, `phase_project.paths(projectRoot)` per i dati e `phase_project.runtime(projectRoot)` per entrambi. Le tre interfacce HTML vengono caricate da una cache utente scrivibile (`prefdir/PHASE/ui`), dove risiedono anche cache delle mappe e log live. I file sorgente restano nell'installazione.

## Importazione legacy

`PHASE_ImportLegacy(oldRoot,newRoot)` cerca `output_###` del modulo Model e `ASC_*`/`DSC_*`/`DES_*` con una cartella `EXPORT`. Copia report XLSX, figure, file GIS/MAT finali ed esportazioni StaMPS. Non tocca l'origine. Scrive un catalogo con percorso originale, destinazione relativa, dimensione e stato della copia. Le immagini SAR, le cartelle SNAP e gli intermedi StaMPS restano nel workspace originale.

`PHASE_ImportLegacy(oldRoot,newRoot,'reference')` crea solo il catalogo. In questa modalità i risultati dipendono dalla posizione del vecchio workspace. Il catalogo permette di rilevare i riferimenti mancanti; è possibile rieseguire un'importazione in copia verso un altro progetto se si vuole un archivio portabile.

L'importazione non legge né riscrive i valori scientifici e non converte i vecchi file `.mat`: ne preserva i byte. Un nuovo modulo Model potrà leggere gli export CSV/XLSX importati. Un preprocessing o una sessione StaMPS legacy avviata a metà va completata nel vecchio ambiente o ricominciata nel nuovo progetto.

## Uso immediato

```matlab
PHASE_NewProject('/percorso/Progetto', 'Nome progetto')
PHASE_Hub('/percorso/Progetto')
PHASE_ImportLegacy('/percorso/vecchio/PHASE', '/percorso/NuovoProgetto')
phase_project.runtime('/percorso/NuovoProgetto')
phase_project.importedResults('/percorso/NuovoProgetto')
phase_project.selfTest()
phase_preprocessing_beta.projectPathsSelfTest()
PHASE_Preprocessing_beta('/percorso/NuovoProgetto')
phase_stamps_beta.projectSelfTest()
```

## Gate prima di usare il nuovo layout per il processing

1. Preprocessing: provare l'intera pipeline Sentinel-1 e COSMO-SkyMed su stack reali, verificando i file prodotti e la ripresa da step intermedi. La separazione dei percorsi è implementata, ma i test finora non eseguono SNAP.
2. StaMPS: validare un dataset reale sotto `processing/stamps`, inclusi i path scritti in `input_StaMPS.mat`, i file generati da snap2stamps, TRAIN e la ripresa di una sessione.
3. Testare un progetto nuovo e uno importato su Windows e Linux; per macOS Apple Silicon validare la pipeline scientifica e i binari esterni prima di dichiarare supporto completo.

`PHASE_Model(projectRoot)` legge e salva `input_model.mat` nella radice del progetto e crea `output_###` in `results/model`; salva anche una copia del log in `logs`. Le risorse di calcolo, come geoSplinter, restano nell'installazione. Il launcher senza argomenti conserva il comportamento PHASE 6.

`PHASE_StaMPS(projectRoot)` apre il dataset presente in `processing/stamps` (o permette di sceglierlo se sono più di uno). Il runtime risolve gli export SNAP in `processing/preprocessing/INSAR_<data>` e, a elaborazione conclusa, copia i prodotti finali da `EXPORT` a `results/exports/<dataset>`, preservando gli originali usati da StaMPS. Il progetto salva una copia del log in `logs`. Questo adattatore è stato provato con dati sintetici, non con una pipeline PSI completa.

`PHASE_Preprocessing_beta(projectRoot)` apre un progetto PHASE 7. Le immagini e `input_preprocessing.mat` risiedono in `processing/preprocessing`; i dataset creati per StaMPS in `processing/stamps`. Il motore scrive le configurazioni SNAP e gli script wrapper nel progetto, invocando gli script Python e i graph template dall'installazione. I file ASF del downloader sono conservati nella cartella `downloadasf` del progetto, mentre il codice del downloader rimane nell'installazione. Il launcher senza argomenti conserva il layout PHASE 6.

`PHASE_Hub()` apre una singola finestra dalla quale scegliere o creare il progetto. `PHASE_Hub(projectRoot)` apre direttamente il progetto indicato. La navigazione persistente Project, Preprocessing, StaMPS PSI e Displacement Modeling monta le rispettive interfacce nella stessa finestra; i moduli vengono caricati su richiesta. La sezione StaMPS PSI consente di scegliere tra i dataset trovati in `processing/stamps`; prima che il Preprocessing ne abbia prodotto uno, resta disponibile con un messaggio esplicativo. Il cambio progetto è impedito mentre un'elaborazione o un download è in corso. I launcher beta standalone restano funzionanti per confronto e diagnosi.

Questa è ancora una modalità di sviluppo: sono stati verificati apertura UI, percorsi e configurazioni sintetiche, non una pipeline SNAP/StaMPS completa su immagini reali né Windows/Linux. I file finali scientifici non sono da considerare validati solo perché la UI si avvia. Prima di usarla in produzione servono le prove indicate sopra.
