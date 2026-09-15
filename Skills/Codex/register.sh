#!/usr/bin/env bash
# Register the bundled skill and allow its external personal state directory.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec "${PYTHON:-python3}" "$script_dir/scripts/register.py" "$@"
