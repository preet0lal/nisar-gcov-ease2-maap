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
  --bbox 'W S E N'            Optional geographic clip (omit for full frames)
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
keep_intermediate="no"

while (($#)); do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --input_file|--input_dir|--ease2_res|--backscatter|--frequency|--pol_pair|--mode|--filter_window|--mosaic|--mosaic_across_modes|--use_exception_mask|--lia|--min_coverage|--noise_floor|--workers|--tile_workers|--tile_size_km|--native_block_px|--gap_threshold|--lia_dem_scale|--dem_dir|--nisar_dem_vrt|--dem_cache_dir|--earthdata_cookie|--ease2_mat|--bbox|--keep_intermediate)
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
  if [[ ! -f "$source_file" ]]; then
    echo "ERROR: MAAP file input was not staged/found: $input_file" >&2
    echo "Expected a downloaded .h5/.hdf5 file in input/, inputs/ or /inputs." >&2
    exit 2
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
  if [[ -n "$bbox" ]]; then
    bbox_clean="${bbox//,/ }"
    read -r -a bbox_values <<< "$bbox_clean"
    ((${#bbox_values[@]} == 4)) || { echo "ERROR: --bbox requires W S E N" >&2; exit 2; }
    args+=(--bbox "${bbox_values[@]}")
  fi
  echo "[MAAP] Starting NISAR GCOV EASE2; filter=$filter; resolution=${ease2_res}m; power=$backscatter"
  echo "[MAAP] Input: $input_dir"
  echo "[MAAP] Output: $target"
  # Deliberately foreground: MAAP must wait for the worker's process to finish.
  "${python_cmd[@]}" "$processor" "${args[@]}"
done

echo "[MAAP] Finished successfully; results are under $output_dir"
