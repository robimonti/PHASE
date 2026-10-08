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

python_default=$(command -v python3 || true)
python_path=$(ask_path 'Percorso Python 3.10+ con venv:' "$python_default")
matlab_default=$(command -v matlab || true)
matlab_path=$(ask_path 'Percorso eseguibile MATLAB (bin/matlab):' "$matlab_default")
gpt_default=/opt/esa-snap/bin/gpt
if [ ! -x "$gpt_default" ]; then gpt_default=/opt/snap/bin/gpt; fi
gpt_path=$(ask_path 'Percorso eseguibile ESA SNAP gpt:' "$gpt_default")
stamps_path=$(ask_path 'Cartella StaMPS già preparata (facoltativa):' '')
train_path=$(ask_path 'Cartella TRAIN già preparata (facoltativa):' '')

set -- "$python_path" "$bundle_dir/install-phase-unix.py" \
    --source "$bundle_dir/engine" --python "$python_path" \
    --matlab "$matlab_path" --gpt "$gpt_path"
if [ -n "$stamps_path" ]; then set -- "$@" --stamps "$stamps_path"; fi
if [ -n "$train_path" ]; then set -- "$@" --train "$train_path"; fi

log_file=$(mktemp)
trap 'rm -f -- "$log_file"' EXIT HUP INT TERM
"$@" >"$log_file" 2>&1 &
install_pid=$!
if [ "$dialog" = zenity ]; then
    (while kill -0 "$install_pid" 2>/dev/null; do printf '50\n'; sleep 1; done) |
        zenity --progress --pulsate --auto-close --no-cancel \
            --title='PHASE Installer' --text='Installazione in corso...' || true
fi
if wait "$install_pid"; then
    show_info 'PHASE è installato. Avvialo dal menu applicazioni.'
else
    show_error "Installazione non riuscita. Dettagli:\n$(tail -n 25 "$log_file")"
    exit 1
fi
