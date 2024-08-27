#!/usr/bin/env bash
source "$HOME/miniforge3/bin/activate" "cpdb"

python "${1}" \
    --input_dir $2


