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

# run_maap.sh uses the BASE image Python (python3) for maap_s3.py, which
# downloads the GCOV granule and gets Earthdata S3 credentials for the LIA DEM.
# The base image may lack s3fs (and possibly maap-py): install what is missing.
python3 - <<'PY'
import importlib, importlib.util, subprocess, sys
packages = {"maap": "maap-py", "s3fs": "s3fs"}
missing = [pkg for mod, pkg in packages.items() if importlib.util.find_spec(mod) is None]
if missing:
    print("Installing into base Python", sys.executable, ":", missing, flush=True)
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--no-cache-dir", *missing])
for mod in packages:
    importlib.import_module(mod)
from maap.maap import MAAP  # noqa: F401
print("Base Python ready for maap_s3.py:", sys.executable, flush=True)
PY
