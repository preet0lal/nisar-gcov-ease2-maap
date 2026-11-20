#!/usr/bin/env bash
# MAAP DPS entry point for the standalone NISAR GCOV -> EASE2 processor.
# Unlike the original nohup launcher, MAAP must execute in the FOREGROUND.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
output_dir="${PWD}/output"
mkdir -p "$output_dir"

usage() {
  cat <<'HELP'
Usage: run_maap.sh [--input_file GCOV.h5 | --input_dir DIR] [named options]
   or: run_maap.sh INPUT_FILE [EASE2_RES MODE FILTER_WINDOW LIA FREQUENCY \
                   BACKSCATTER CLIP_TO_AOI BBOX_W BBOX_S BBOX_E BBOX_N]

Positional form (used by MAAP DPS, which passes registered inputs as bare
values in registration order). Empty, "none" or "-" keeps the default.
CLIP_TO_AOI=yes crops outputs to BBOX_W BBOX_S BBOX_E BBOX_N (EPSG:4326);
CLIP_TO_AOI=no keeps full frames and ignores the bbox values.

MAAP DPS automatically downloads a registered File Input into the job's input
folder. --input_file accepts that staged path, the source URL (resolved by
basename in input/, inputs/ or /inputs), or an absolute local path. If omitted,
all GCOV .h5/.hdf5 in --input_dir or a known input folder are discovered.

Options:
  --input_file PATH_OR_URI    One MAAP-staged NISAR GCOV HDF5 granule
  --input_dir PATH            Directory of GCOV granules (for local/batch use)
  --ease2_res 50|100          Default 100 (M0p1)
  --backscatter gamma0|sigma0 Default gamma0, LINEAR output power
  --frequency frequencyA|frequencyB     Default frequencyA
  --pol_pair auto|HH_HV|VV_VH|RH_RV      Default auto
  --mode nofilt|enhanced_lee|gamma_map|all  Default nofilt
  --filter_window N           Default 7
  --mosaic yes|no             Default yes
  --mosaic_across_modes yes|no           Default no
  --use_exception_mask yes|no           Default yes (CRID-aware)
  --lia yes|no                Default no
  --min_coverage FRACTION     Default 0.90
  --noise_floor POWER         Default 0
  --workers N                 Default 1
  --tile_workers N            Default 4
  --tile_size_km N            Default 50
  --native_block_px N         Default 1024
  --gap_threshold KM          Default 50
  --lia_dem_scale METERS      Default 30
  --dem_dir DIR               Optional local DEM tile folder, for LIA
  --nisar_dem_vrt PATH_OR_URL Optional DEM VRT, for LIA
  --dem_cache_dir DIR         Optional writable cache, for LIA
  --earthdata_cookie PATH     Optional authenticated cookie file, for LIA
  --ease2_mat PATH            Optional grid-coordinate validation file
  --clip_to_aoi yes|no        Crop outputs to the AOI bbox (default: yes if a bbox is given)
  --bbox 'W S E N'            AOI bbox in EPSG:4326 (or --bbox_w/--bbox_s/--bbox_e/--bbox_n)
  --keep_intermediate yes|no   Default no
  --help

Outputs are ALWAYS saved in the MAAP DPS 'output/' directory.
HELP
}

input_file=""
input_dir=""
ease2_res="100"
backscatter="gamma0"
frequency="frequencyA"
pol_pair="auto"
mode="nofilt"
filter_window="7"
mosaic="yes"
mosaic_across_modes="no"
use_exception_mask="yes"
lia="no"
min_coverage="0.90"
noise_floor="0"
workers="1"
tile_workers="4"
tile_size_km="50"
native_block_px="1024"
gap_threshold="50"
lia_dem_scale="30"
dem_dir=""
nisar_dem_vrt=""
dem_cache_dir=""
earthdata_cookie=""
ease2_mat=""
bbox=""
clip_to_aoi=""
bbox_w=""
bbox_s=""
bbox_e=""
bbox_n=""
keep_intermediate="no"

# Positional (MAAP DPS) form: first argument is not a --flag.
if (($#)) && [[ "$1" != -* ]]; then
  positional_names=(input_file ease2_res mode filter_window lia frequency
                    backscatter clip_to_aoi bbox_w bbox_s bbox_e bbox_n)
  (($# <= ${#positional_names[@]})) || { echo "ERROR: too many positional arguments ($#)" >&2; exit 2; }
  index=0
  for value in "$@"; do
    case "${value,,}" in ""|none|-|null) ;; *) printf -v "${positional_names[$index]}" '%s' "$value" ;; esac
    index=$((index + 1))
  done
  set --
fi

while (($#)); do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --input_file|--input_dir|--ease2_res|--backscatter|--frequency|--pol_pair|--mode|--filter_window|--mosaic|--mosaic_across_modes|--use_exception_mask|--lia|--min_coverage|--noise_floor|--workers|--tile_workers|--tile_size_km|--native_block_px|--gap_threshold|--lia_dem_scale|--dem_dir|--nisar_dem_vrt|--dem_cache_dir|--earthdata_cookie|--ease2_mat|--bbox|--clip_to_aoi|--bbox_w|--bbox_s|--bbox_e|--bbox_n|--keep_intermediate)
      (($# >= 2)) || { echo "Missing value for $1" >&2; exit 2; }
      key="${1#--}"
      value="$2"
      printf -v "$key" '%s' "$value"
      shift 2 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

# An explicitly supplied file means exactly ONE granule, even if the MAAP job
# input/ folder holds other products. Symlink rather than copy multi-GB HDF5.
base_python="${MAAP_BASE_PYTHON:-python3}"   # MAAP base image Python (maap-py, s3fs)
stage_dir=""
cleanup() { if [[ -n "$stage_dir" && -d "$stage_dir" ]]; then rm -rf -- "$stage_dir"; fi; }
trap cleanup EXIT

if [[ -n "$input_file" ]]; then
  source_file=""
  if [[ -f "$input_file" ]]; then
    source_file="$input_file"
  else
    # DPS file input may pass an s3/https URL even though the file is staged.
    filename="${input_file%%\?*}"
    filename="${filename##*/}"
    for directory in "$PWD/input" "$PWD/inputs" /inputs "$PWD"; do
      [[ -d "$directory" ]] || continue
      candidate="$(find "$directory" -maxdepth 4 -type f -name "$filename" -print -quit 2>/dev/null || true)"
      if [[ -n "$candidate" ]]; then source_file="$candidate"; break; fi
    done
  fi
  if [[ ! -f "$source_file" && "$input_file" == s3://* ]]; then
    # Not staged by DPS: fetch it with temporary Earthdata S3 credentials.
    source_file="$PWD/input/$filename"
    if [[ "$input_file" == s3://maap-ops-workspace/* ]]; then
      # A CSV manifest in the MAAP private workspace bucket may not be staged
      # automatically by older DPS builds; use workspace credentials if needed.
      "$base_python" "${script_dir}/maap_workspace_s3.py" get "$input_file" "$source_file"
    else
      "$base_python" "${script_dir}/maap_s3.py" get "$input_file" "$source_file"
    fi
  fi
  if [[ ! -f "$source_file" ]]; then
    echo "ERROR: MAAP file input was not staged/found: $input_file" >&2
    echo "Expected a downloaded .h5/.hdf5 file in input/, inputs/ or /inputs." >&2
    exit 2
  fi
  # Manifest CSV => one DPS job sequentially processes EVERY listed full GCOV frame.
  # Individual .h5/.hdf5 inputs retain the previous one-scene behavior.
  if [[ "${source_file,,}" == *.csv ]]; then
    echo "[MAAP] Manifest file detected: processing all GCOV scenes in ONE MAAP DPS job"
    exec "$base_python" "${script_dir}/maap_batch_manifest.py" \
      --manifest "$source_file" \
      --runner "${script_dir}/run_maap.sh" \
      --downloader "${script_dir}/maap_s3.py" \
      --output-dir "$output_dir" \
      --base-python "$base_python" \
      --ease2-res "$ease2_res" \
      --backscatter "$backscatter" \
      --mode "$mode" \
      --filter-window "$filter_window" \
      --lia "$lia" \
      --frequency "$frequency" \
      --pol-pair "$pol_pair" \
      --workers "$workers" \
      --tile-workers "$tile_workers"
  fi
  case "${source_file,,}" in *.h5|*.hdf5) ;; *) echo "ERROR: Expected a .h5 or .hdf5: $source_file" >&2; exit 2 ;; esac
  stage_dir="$(mktemp -d "${TMPDIR:-/tmp}/maap_gcov_ease2_XXXXXXXX")"
  ln -s "$(realpath "$source_file")" "$stage_dir/$(basename "$source_file")"
  input_dir="$stage_dir"
elif [[ -n "$input_dir" ]]; then
  [[ -d "$input_dir" ]] || { echo "ERROR: --input_dir not found: $input_dir" >&2; exit 2; }
else
  for directory in "$PWD/input" "$PWD/inputs" /inputs; do
    [[ -d "$directory" ]] || continue
    if find "$directory" -type f \( -iname '*.h5' -o -iname '*.hdf5' \) -print -quit | grep -q .; then
      input_dir="$directory"
      break
    fi
  done
  [[ -n "$input_dir" ]] || { echo "ERROR: No GCOV inputs found. Supply --input_file or --input_dir." >&2; exit 2; }
fi

case "$mode" in
  nofilt) modes=(none) ;;
  enhanced_lee) modes=(enhanced_lee) ;;
  gamma_map) modes=(gamma_map) ;;
  all) modes=(none enhanced_lee gamma_map) ;;
  *) echo "ERROR: Invalid --mode: $mode" >&2; exit 2 ;;
esac

if [[ -n "${MAAP_PYTHON_EXEC:-}" ]]; then
  python_cmd=("$MAAP_PYTHON_EXEC")
elif command -v conda >/dev/null 2>&1; then
  python_cmd=(conda run --live-stream --name nisar_gcov_ease2 python)
else
  echo "ERROR: conda not available; run build-env.sh during MAAP registration." >&2
  exit 1
fi

# ---- AOI clipping ----------------------------------------------------------
if [[ -z "$bbox" && -n "$bbox_w$bbox_s$bbox_e$bbox_n" ]]; then
  bbox="$bbox_w $bbox_s $bbox_e $bbox_n"
fi
[[ -n "$clip_to_aoi" ]] || { [[ -n "$bbox" ]] && clip_to_aoi="yes" || clip_to_aoi="no"; }
case "$clip_to_aoi" in yes|no) ;; *) echo "ERROR: --clip_to_aoi must be yes or no" >&2; exit 2 ;; esac
bbox_values=()
if [[ "$clip_to_aoi" == "yes" ]]; then
  read -r -a bbox_values <<< "${bbox//,/ }"
  ((${#bbox_values[@]} == 4)) || { echo "ERROR: clip_to_aoi=yes requires a bbox W S E N" >&2; exit 2; }
  python3 - "${bbox_values[@]}" <<'PY' || exit 2
import sys
try:
    w, s, e, n = map(float, sys.argv[1:])
except ValueError:
    sys.exit("ERROR: bbox values must be numbers: " + " ".join(sys.argv[1:]))
if not (-180 <= w < e <= 180 and -90 <= s < n <= 90):
    sys.exit(f"ERROR: invalid bbox W={w} S={s} E={e} N={n} (need W<E, S<N)")
PY
  echo "[MAAP] AOI clip: outputs cropped to ${bbox_values[*]} (EPSG:4326)"
else
  echo "[MAAP] AOI clip: no (full frames)"
fi

# ---- LIA DEM access on DPS -------------------------------------------------
# The NISAR DEM needs Earthdata access; on a DPS worker use temporary S3
# credentials and the in-region /vsis3 path unless a DEM source was given.
if [[ "$lia" == "yes" && -z "$dem_dir" && -z "$nisar_dem_vrt" ]]; then
  if creds="$("$base_python" "${script_dir}/maap_s3.py" env 2>/dev/null)"; then
    eval "$creds"; unset creds
    nisar_dem_vrt="/vsis3/sds-n-cumulus-prod-nisar-products/DEM/v1.2/EPSG4326/EPSG4326.vrt"
    echo "[MAAP] LIA DEM: $nisar_dem_vrt (temporary Earthdata S3 credentials)"
  else
    echo "[MAAP] WARNING: could not get Earthdata S3 credentials; LIA uses the processor's default DEM URL" >&2
  fi
fi

processor="${MAAP_PROCESSOR_PATH:-${script_dir}/NISAR_GCOV_TO_EASE2.py}"
[[ -f "$processor" ]] || { echo "ERROR: Python processor missing: $processor" >&2; exit 1; }

for filter in "${modes[@]}"; do
  target="$output_dir"
  if ((${#modes[@]} > 1)); then
    subdir="$filter"; [[ "$subdir" == "none" ]] && subdir="nofilt"
    target="$output_dir/$subdir"
  fi
  mkdir -p "$target"
  args=(
    --input_dir "$input_dir" --output_dir "$target"
    --ease2_res "$ease2_res" --backscatter "$backscatter"
    --frequency "$frequency" --pol_pair "$pol_pair"
    --filter "$filter" --filter_window "$filter_window"
    --mosaic "$mosaic" --mosaic_across_modes "$mosaic_across_modes"
    --use_exception_mask "$use_exception_mask" --lia "$lia"
    --min_coverage "$min_coverage" --noise_floor "$noise_floor"
    --workers "$workers" --tile_workers "$tile_workers"
    --tile_size_km "$tile_size_km" --native_block_px "$native_block_px"
    --gap_threshold "$gap_threshold" --lia_dem_scale "$lia_dem_scale"
    --keep_intermediate "$keep_intermediate"
  )
  [[ -z "$dem_dir" ]] || args+=(--dem_dir "$dem_dir")
  [[ -z "$nisar_dem_vrt" ]] || args+=(--nisar_dem_vrt "$nisar_dem_vrt")
  [[ -z "$dem_cache_dir" ]] || args+=(--dem_cache_dir "$dem_cache_dir")
  [[ -z "$earthdata_cookie" ]] || args+=(--earthdata_cookie "$earthdata_cookie")
  [[ -z "$ease2_mat" ]] || args+=(--ease2_mat "$ease2_mat")
  if [[ "$clip_to_aoi" == "yes" ]]; then
    args+=(--bbox "${bbox_values[@]}")
  fi
  echo "[MAAP] Starting NISAR GCOV EASE2; filter=$filter (window $filter_window); resolution=${ease2_res}m; power=$backscatter; frequency=$frequency; lia=$lia; clip_to_aoi=$clip_to_aoi"
  echo "[MAAP] Input: $input_dir"
  echo "[MAAP] Output: $target"
  # Deliberately foreground: MAAP must wait for the worker's process to finish.
  "${python_cmd[@]}" "$processor" "${args[@]}"
done

echo "[MAAP] Finished successfully; results are under $output_dir"
