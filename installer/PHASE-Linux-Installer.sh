#!/bin/sh
set -eu

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
bundle_dir=$(CDPATH= cd -- "$script_dir/../share/phase" && pwd)

if command -v zenity >/dev/null 2>&1; then
    dialog=zenity
elif command -v kdialog >/dev/null 2>&1; then
    dialog=kdialog
else
    printf '%s\n' 'PHASE installer needs zenity or kdialog for its graphical setup.' >&2
    exit 1
fi

show_error() {
    if [ "$dialog" = zenity ]; then
        zenity --error --title='PHASE Installer' --text="$1"
    else
        kdialog --error "$1" --title 'PHASE Installer'
    fi
}

show_info() {
    if [ "$dialog" = zenity ]; then
        zenity --info --title='PHASE Installer' --text="$1"
    else
        kdialog --msgbox "$1" --title 'PHASE Installer'
    fi
}

ask_path() {
    if [ "$dialog" = zenity ]; then
        zenity --entry --title='PHASE Installer' --text="$1" --entry-text="$2"
    else
        kdialog --inputbox "$1" "$2" --title 'PHASE Installer'
    fi
}

show_info 'PHASE is installed once for all projects. Create or open project folders anywhere after launching the app. MATLAB, SNAP, Python and Linux build tools are required.'
python_default=$(command -v python3 || true)
python_path=$(ask_path 'Python 3.10+ executable (with venv):' "$python_default")
matlab_default=$(command -v matlab || true)
matlab_path=$(ask_path 'MATLAB executable (bin/matlab):' "$matlab_default")
gpt_default=/opt/esa-snap/bin/gpt
if [ ! -x "$gpt_default" ]; then gpt_default=/opt/snap/bin/gpt; fi
gpt_path=$(ask_path 'ESA SNAP gpt executable:' "$gpt_default")

runtime_tmp=$(mktemp -d "${TMPDIR:-/tmp}/phase-linux-installer.XXXXXXXX")
log_file="$runtime_tmp/install.log"
cleanup() {
    case "$runtime_tmp" in
        "${TMPDIR:-/tmp}"/phase-linux-installer.*)
            if [ -d "$runtime_tmp" ]; then rm -r -- "$runtime_tmp"; fi ;;
    esac
}
trap cleanup EXIT HUP INT TERM

"$python_path" "$bundle_dir/prepare-linux-runtime.py" \
    --output "$runtime_tmp/runtime" >"$log_file" 2>&1 &
build_pid=$!
if [ "$dialog" = zenity ]; then
    (while kill -0 "$build_pid" 2>/dev/null; do printf '50\n'; sleep 1; done) |
        zenity --progress --pulsate --auto-close --no-cancel \
            --title='PHASE Installer' --text='Preparing StaMPS and TRAIN…' || true
fi
if ! wait "$build_pid"; then
    show_error "StaMPS/TRAIN preparation failed. Details:\n$(tail -n 25 "$log_file")"
    exit 1
fi

set -- "$python_path" "$bundle_dir/install-phase-unix.py" \
    --source "$bundle_dir/engine" --python "$python_path" \
    --matlab "$matlab_path" --gpt "$gpt_path" \
    --stamps "$runtime_tmp/runtime/StaMPS" \
    --train "$runtime_tmp/runtime/TRAIN"

"$@" >"$log_file" 2>&1 &
install_pid=$!
if [ "$dialog" = zenity ]; then
    (while kill -0 "$install_pid" 2>/dev/null; do printf '50\n'; sleep 1; done) |
        zenity --progress --pulsate --auto-close --no-cancel \
            --title='PHASE Installer' --text='Installing PHASE…' || true
fi
if wait "$install_pid"; then
    if [ "$dialog" = zenity ]; then
        if zenity --question --title='PHASE Installer' \
            --text='Installation complete. Launch PHASE now?'; then
            "$HOME/.local/share/PHASE/launch-phase.sh" >/dev/null 2>&1 &
        fi
    elif kdialog --yesno 'Installation complete. Launch PHASE now?' \
            --title 'PHASE Installer'; then
        "$HOME/.local/share/PHASE/launch-phase.sh" >/dev/null 2>&1 &
    fi
else
    show_error "Installation failed. Details:\n$(tail -n 25 "$log_file")"
    exit 1
fi
