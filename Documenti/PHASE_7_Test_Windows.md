# PHASE 7 — prova end-to-end su Windows

Questa prova va eseguita dopo che il branch PHASE 7 è disponibile su GitHub e
l'installer Windows è stato compilato dal medesimo codice. L'EXE pubblico 6.1.4
non contiene l'hub. Usare una copia di uno stack SAR noto, con spazio disco
sufficiente per SNAP e StaMPS; conservare intatto il dataset originale.

Per questa anteprima, aprire su GitHub la sezione **Actions**, scegliere
**Build PHASE 7 Windows test installer**, aprire l'ultima esecuzione riuscita
del branch `codex/phase-stamps-beta` e scaricare l'artefatto
`phase7-windows-test-installer`. GitHub lo consegna in uno ZIP che contiene un
solo `install-phase7-test.exe`: estrarlo e avviare l'EXE. Il branch di test è
incorporato nel wizard, quindi non servono parametri da terminale. Il pacchetto
non è una release pubblica e può mostrare un avviso SmartScreen perché non è
firmato.

Per aggiornare una precedente installazione di test, chiudere PHASE/MATLAB e
rieseguire l'EXE scegliendo la stessa destinazione. L'installer aggiorna il
codice nell'installazione condivisa senza rimuovere le cartelle dei progetti;
rinomina il proprio collegamento desktop da `PHASE 7` a `PHASE`.

## Preparazione

1. Verificare MATLAB, licenza, ESA SNAP 13 e accesso a GitHub.
2. Lanciare il nuovo installer Windows. Scegliere una destinazione utente
   (predefinita `%LOCALAPPDATA%\Programs\PHASE`).
3. Verificare `PHASE.lnk` sul desktop con il logo PHASE,
   `PHASE\PHASE.lnk`, `PHASE\install.json`,
   `PHASE\engine\PHASE_Hub.m`, `PHASE\engine\StaMPS` e `PHASE\engine\TRAIN`.
4. Avviare PHASE dal collegamento, non da una vecchia scorciatoia 6.1.
   Se l'installazione fallisce, usare **Apri log** nel wizard: il file si trova
   in `%TEMP%\phase-installer-*.log`.

## Progetto e pipeline

1. Creare un progetto nuovo in una cartella diversa da quella di installazione.
   Verificare `phase-project.json`, `input`, `processing`, `results` e `logs`.
2. Nella scheda Preprocessing impostare AOI e uno stack Sentinel-1 di prova.
   Eseguire tutte le fasi SNAP. Controllare input e intermedi in
   `processing\preprocessing`; annotare eventuali errori del monitor.
3. Aprire la scheda StaMPS e selezionare il dataset creato in
   `processing\stamps\ASC_*` o `DSC_*`. Eseguire la sequenza PSI completa,
   compresi i binari nativi. Se si usa TRAIN, verificare separatamente la
   correzione atmosferica. Controllare `results\exports` e i log del progetto.
4. Nella scheda Model caricare gli export, avviare un'elaborazione completa e
   controllare report, figure e file GIS in `results\model` e nelle cartelle
   `results` previste dal modulo.
5. Chiudere PHASE e MATLAB, riaprire il progetto dal collegamento e verificare
   che configurazioni, dataset ed export siano ancora disponibili.

## Aggiornamento e risultato della prova

Il pulsante **Cerca update** può essere provato solo quando esiste una release
stabile v7 con asset `phase7-engine.zip`. Dopo il download, chiudere MATLAB e
riavviare dal collegamento. Controllare la versione in `install.json`, il backup
in `backups` e la persistenza del progetto e dei runtime StaMPS/TRAIN.

Per ogni eventuale errore registrare la fase, il messaggio completo, la
versione PHASE/SNAP/MATLAB, il log pertinente e il percorso del file atteso.
Non considerare superata la prova solo perché l'hub e le schede si aprono:
occorrono i prodotti finali scientifici e un controllo di plausibilità su uno
stack già conosciuto.
