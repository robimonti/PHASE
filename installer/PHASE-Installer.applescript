use scripting additions

on run
    set bundleRoot to POSIX path of (path to me)
    set resourcesRoot to bundleRoot & "Contents/Resources/"
    set installerScript to resourcesRoot & "install-phase-unix.py"
    set engineSource to resourcesRoot & "engine"

    try
        set pythonPath to do shell script ("sh " & quoted form of (resourcesRoot & "find-python-macos.sh"))
    on error
        display dialog "Per installare PHASE serve Python 3.10 o più recente con venv. Se è installato, seleziona il suo eseguibile; altrimenti installalo e riapri questo installer." buttons {"Seleziona", "Annulla"} default button "Seleziona" with icon caution
        set pythonPath to POSIX path of (choose file with prompt "Eseguibile Python 3.10+")
    end try

    set matlabPath to do shell script "ls -1d /Applications/MATLAB_R*.app/bin/matlab 2>/dev/null | tail -n 1"
    if matlabPath is "" then
        display dialog "Seleziona l'eseguibile matlab dentro MATLAB_R*.app/bin." buttons {"Seleziona", "Annulla"} default button "Seleziona"
        set matlabPath to POSIX path of (choose file with prompt "Eseguibile MATLAB")
    end if

    set snapPath to "/Applications/esa-snap/bin/gpt"
    try
        do shell script "test -x " & quoted form of snapPath
    on error
        display dialog "Seleziona l'eseguibile gpt di ESA SNAP." buttons {"Seleziona", "Annulla"} default button "Seleziona"
        set snapPath to POSIX path of (choose file with prompt "Eseguibile SNAP gpt")
    end try

    set answer to display dialog "PHASE sarà installato nella tua cartella utente e apparirà in Applicazioni. MATLAB e SNAP resteranno nelle loro installazioni attuali. Continuare?" buttons {"Annulla", "Installa PHASE"} default button "Installa PHASE" with icon note
    if button returned of answer is not "Installa PHASE" then return

    display notification "L'installazione può richiedere alcuni minuti." with title "PHASE"
    set commandText to quoted form of pythonPath & " " & quoted form of installerScript & " --source " & quoted form of engineSource & " --python " & quoted form of pythonPath & " --matlab " & quoted form of matlabPath & " --gpt " & quoted form of snapPath
    try
        with timeout of 3600 seconds
            set installResult to do shell script commandText & " 2>&1"
        end timeout
        display dialog "PHASE è installato. Apri PHASE.app dalla cartella Applicazioni del tuo utente.\n\n" & installResult buttons {"OK"} default button "OK" with icon note
    on error messageText
        display dialog "Installazione non riuscita:\n\n" & messageText buttons {"OK"} default button "OK" with icon caution
    end try
end run
