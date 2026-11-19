#!/usr/bin/env bash
# MAAP DPS build command. Create an isolated GDAL/Rasterio-compatible environment.
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda was not found in the selected MAAP base container" >&2
  exit 1
fi
if conda env list | awk '{print $1}' | grep -qx 'nisar_gcov_ease2'; then
  conda env update --name nisar_gcov_ease2 --file "${script_dir}/environment.yml" --prune
else
  conda env create --file "${script_dir}/environment.yml"
fi
conda run --live-stream --name nisar_gcov_ease2 python -c \
  'import h5py,numpy,rasterio,scipy,matplotlib; print("MAAP EASE2 environment ready; rasterio",rasterio.__version__)'
