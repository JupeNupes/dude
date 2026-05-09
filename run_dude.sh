#!/bin/bash
# Launch script for dude.py
# Required because GTK libraries are installed via Homebrew at /opt/homebrew/lib
# and macOS does not search there by default.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export DYLD_LIBRARY_PATH="/opt/homebrew/lib:$DYLD_LIBRARY_PATH"

source "$SCRIPT_DIR/.venv/bin/activate"
python "$SCRIPT_DIR/dude.py"
