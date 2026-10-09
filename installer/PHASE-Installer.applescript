use scripting additions

on run
    set bundleRoot to POSIX path of (path to me)
    set resourcesRoot to bundleRoot & "Contents/Resources/"
    set installerScript to resourcesRoot & "install-phase-unix.py"
    set engineSource to resourcesRoot & "engine"

    try
        set pythonPath to do shell script ("sh " & quoted form of (resourcesRoot & "find-python-macos.sh"))
    on error
        display dialog "PHASE needs Python 3.10 or newer with venv. Select its executable, or install Python and reopen this installer." buttons {"Select", "Cancel"} default button "Select" with icon caution
        set pythonPath to POSIX path of (choose file with prompt "Python 3.10+ executable")
    end try

    set matlabPath to do shell script "ls -1d /Applications/MATLAB_R*.app/bin/matlab 2>/dev/null | tail -n 1"
    if matlabPath is "" then
        display dialog "Select the matlab executable inside MATLAB_R*.app/bin." buttons {"Select", "Cancel"} default button "Select"
        set matlabPath to POSIX path of (choose file with prompt "MATLAB executable")
    end if

    set snapPath to "/Applications/esa-snap/bin/gpt"
    try
        do shell script "test -x " & quoted form of snapPath
    on error
        display dialog "Select the ESA SNAP gpt executable." buttons {"Select", "Cancel"} default button "Select"
        set snapPath to POSIX path of (choose file with prompt "SNAP gpt executable")
    end try

    set runtimeArgs to ""
    try
        do shell script "test -f " & quoted form of (resourcesRoot & "StaMPS/matlab/stamps.m") & " -a -f " & quoted form of (resourcesRoot & "TRAIN/matlab/aps_linear.m")
        set runtimeArgs to " --stamps " & quoted form of (resourcesRoot & "StaMPS") & " --train " & quoted form of (resourcesRoot & "TRAIN")
        set psiMessage to "This installer includes the Apple Silicon StaMPS and TRAIN runtime."
    on error
        set psiMessage to "This installer does not include a native PSI runtime. Preprocessing and Modeling can be installed, but StaMPS processing will need a prepared runtime."
    end try

    set answer to display dialog "PHASE is installed once in your user Applications folder. You can create any number of projects in folders anywhere you choose; project data stays separate from the app. MATLAB and SNAP remain in their current locations.\n\n" & psiMessage buttons {"Cancel", "Install PHASE"} default button "Install PHASE" with icon note
    if button returned of answer is not "Install PHASE" then return

    display notification "Installation may take a few minutes." with title "PHASE"
    set commandText to quoted form of pythonPath & " " & quoted form of installerScript & " --source " & quoted form of engineSource & " --python " & quoted form of pythonPath & " --matlab " & quoted form of matlabPath & " --gpt " & quoted form of snapPath & runtimeArgs
    try
        with timeout of 3600 seconds
            set installResult to do shell script commandText & " 2>&1"
        end timeout
        set finishAnswer to display dialog "PHASE is installed. Open PHASE.app from your user Applications folder, or launch it now.\n\n" & installResult buttons {"Close", "Launch PHASE"} default button "Launch PHASE" with icon note
        if button returned of finishAnswer is "Launch PHASE" then
            set appPath to (POSIX path of (path to home folder)) & "Library/Application Support/PHASE/PHASE.app"
            try
                do shell script "open " & quoted form of appPath
            on error launchError
                display dialog "PHASE was installed, but could not launch automatically:\n\n" & launchError buttons {"OK"} default button "OK" with icon caution
            end try
        end if
    on error messageText
        display dialog "Installation failed:\n\n" & messageText buttons {"OK"} default button "OK" with icon caution
    end try
end run
