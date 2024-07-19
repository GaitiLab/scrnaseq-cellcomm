#!/usr/bin/env bash
source "$HOME/miniforge3/bin/activate" "cpdb"

python "${1}/Python/update_cellphonedb.py" \
    --input_dir $2


