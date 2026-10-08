#!/bin/sh
set -eu

for candidate in \
    /Library/Frameworks/Python.framework/Versions/Current/bin/python3 \
    /opt/homebrew/bin/python3 \
    /usr/local/bin/python3 \
    /usr/bin/python3
do
    if [ -x "$candidate" ] && "$candidate" -c \
        'import sys, venv; sys.exit(0 if sys.version_info >= (3, 10) else 1)' \
        >/dev/null 2>&1; then
        printf '%s\n' "$candidate"
        exit 0
    fi
done
exit 1
