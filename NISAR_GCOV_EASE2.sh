#!/usr/bin/env bash
# NISAR_GCOV_EASE2.sh -- launcher for NISAR_GCOV_TO_EASE2.py
# Last updated: 2026-11-19
#
# Change log
#   2026-11-19  First version. EASE-Grid 2.0 only (50 m or 100 m) with
#               gamma0 / sigma0 output. Derived from NISAR_GCOV_TIFF.sh
#               (2026-11-18); that launcher and its processor are unchanged.
#               Default Earthdata cookie is ~/.gdal_edl_cookies (a stale
#               shared /tmp/gdal_cookies.txt broke DEM downloads on gsapp17).
set -euo pipefail

SCRIPT_VERSION="2026-11-19"

usage() {
cat <<'USAGE'
NISAR_GCOV_EASE2.sh  (last updated 2026-11-19)

Usage:
  bash NISAR_GCOV_EASE2.sh \
    --input_dir INPUT \
    --output_dir OUTPUT \
    --res 50|100 \
    --backscatter gamma0|sigma0 \
    [options]

Purpose:
  Convert NISAR L2 GCOV HDF5 files to GeoTIFFs on the global EASE-Grid 2.0
  (EPSG:6933), at 50 m or 100 m, as gamma-naught or sigma-naught.
  Values are LINEAR power; quicklooks show dB.

Grid (--res):
  100  M0p1 : 100.089502334956 m, 347040 cols x 146160 rows (global)
   50  M0p05:  50.044751167478 m, 694080 cols x 292320 rows (global)
  Both nest exactly in the 1 km EASE2 M01 grid: one 1 km cell = 10x10 M0p1
  cells = 20x20 M0p05 cells, and one 100 m cell = exactly 2x2 50 m cells.
  Each output gets a *_EASE2_grid.mat with global zero-based row/col, the
  parent 1 km M01 row/col, and x/y/lon/lat pixel-centre vectors.

Radiometry (--backscatter):
  gamma0  RTC gamma-naught as stored in GCOV (default)
  sigma0  sigma-naught = gamma0 x rtcGammaToSigmaFactor, applied to every
          native pixel BEFORE speckle filtering and aggregation. The factor
          is ~cos(local incidence angle) (median ~0.76), so sigma0 is about
          1.2 dB below gamma0 on flat ground at 40 deg. The product's stored
          convention is checked; nothing is converted twice.

Output names:
  NISAR_L_GCOV_c<cycle>_t<track>_<pass>_f<frames>_<mode>_<pol>_<freq>_
  <date>_<CRID>_EASE2_<50|100>m_<gamma0|sigma0>_<filter>_<term>.tif

Product maturity by CRID (unchanged from NISAR_GCOV_TO_TIFF.py):
  CRID <  5030 provisional: exception mask uint8, valid == 0 (frequencyB
               layer folded 8x8 from the 10 m grid)
  CRID >= 5030 validated: uint16; not fill, low byte 0, both pol validity
               bits set (good HH/HV = 768)

Required:
  --input_dir PATH          Folder containing NISAR GCOV .h5/.hdf5 files
  --output_dir PATH         Output folder

Grid / radiometry:
  --res 50|100              EASE2 resolution (default 100); alias --ease2_res
  --backscatter TYPE        gamma0 | sigma0 (default gamma0);
                            alias --normalization
  --ease2_mat PATH          Optional lat/lon MAT check of every pixel centre,
                            e.g. EASE2Inverse_latlon50m.mat for --res 50
  --bbox "W S E N"          Optional clip box in EPSG:4326

Data options:
  --frequency NAME          frequencyA | frequencyB (default frequencyA).
                            frequencyB is posted at 80 m; at --res 50 it is
                            oversampled (the log warns).
  --pol_pair PAIR           auto | HH_HV | VV_VH | RH_RV (default auto)
  --mode MODE               nofilt | enhanced_lee | gamma_map | all
                            (default nofilt)
  --window N                Speckle-filter window, odd >= 3 (default 7)
  --noise_floor VALUE       Linear-power threshold; 0 disables (default 0)
  --min_coverage F          Minimum valid native area per output pixel
                            (default 0.90)
  --use_exception_mask VAL  yes | no (default yes)
  --lia yes|no              Incidence + local incidence angle GeoTIFFs on the
                            same grid (default no)
  --lia_dem_scale M         Internal DEM scale for LIA (default 30 m)
  --dem_dir PATH            Local NISAR DEM tiles (no download needed)
  --nisar_dem_vrt PATH      DEM VRT override (default NISAR DEM v1.2 on ASF)
  --dem_cache_dir PATH      DEM cache (default <input_dir>/DEM_Cache)
  --earthdata_cookie PATH   GDAL cookie seed (default ~/.gdal_edl_cookies)

Mosaic/performance options:
  --mosaic yes|no           Mosaic adjacent compatible frames (default yes)
  --mosaic_across_modes yes|no
                            Merge 2005 and 4005 frames (default no; 2005 reads
                            ~0.25 dB HH / ~0.18 dB HV lower than 4005)
  --gap_threshold KM        Maximum inter-frame gap (default 50)
  --workers N               Parallel acquisition groups (default 1)
  --tile_workers N          GDAL warp threads per group (default 4)
  --tile_size_km KM         Output processing window size (default 50)
  --native_block_px N       Native HDF5 block size (default 1024)
  --keep_intermediate yes|no

Execution options:
  --script PATH             Python processor path
                            (default: NISAR_GCOV_TO_EASE2.py next to this file)
  --python_exec PATH        Python executable (default python3)
  --extra_args "ARGS"       Extra Python arguments
  -h, --help

Examples:
  # 50 m sigma0, enhanced Lee, with grid check against the 50 m MAT
  bash NISAR_GCOV_EASE2.sh \
    --input_dir /gpfs/data1/vclgp/lalp/NISAR_Data/Kruger \
    --output_dir /gpfs/data1/vclgp/lalp/NISAR_Data/Kruger/EASE2_50m_sigma0 \
    --res 50 --backscatter sigma0 --mode enhanced_lee --window 7 \
    --ease2_mat /path/to/EASE2Inverse_latlon50m.mat

  # 100 m gamma0 with LIA
  bash NISAR_GCOV_EASE2.sh \
    --input_dir /path/to/gcov --output_dir /path/to/out \
    --res 100 --backscatter gamma0 --lia yes
USAGE
}

INPUT_DIR=""
OUTPUT_DIR=""
RES="100"
BACKSCATTER="gamma0"
EASE2_MAT=""
MODE="nofilt"
MOSAIC="yes"
MOSAIC_ACROSS_MODES="no"
WORKERS="1"
TILE_WORKERS="4"
TILE_SIZE_KM="50"
NATIVE_BLOCK_PX="1024"
FREQUENCY="frequencyA"
POL_PAIR="auto"
WINDOW="7"
NOISE_FLOOR="0"
MIN_COVERAGE="0.90"
USE_EXCEPTION_MASK="yes"
LIA="no"
LIA_DEM_SCALE="30"
DEM_DIR=""
NISAR_DEM_VRT=""
DEM_CACHE_DIR=""
EARTHDATA_COOKIE="${HOME}/.gdal_edl_cookies"
GAP_THRESHOLD="50"
BBOX=""
KEEP_INTERMEDIATE="no"
EXTRA_ARGS=""
PYTHON_EXEC="python3"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/NISAR_GCOV_TO_EASE2.py"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input_dir) INPUT_DIR="$2"; shift 2 ;;
    --output_dir) OUTPUT_DIR="$2"; shift 2 ;;
    --res|--ease2_res) RES="$2"; shift 2 ;;
    --backscatter|--normalization) BACKSCATTER="$2"; shift 2 ;;
    --ease2_mat) EASE2_MAT="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --mosaic) MOSAIC="$2"; shift 2 ;;
    --mosaic_across_modes) MOSAIC_ACROSS_MODES="$2"; shift 2 ;;
    --workers) WORKERS="$2"; shift 2 ;;
    --tile_workers) TILE_WORKERS="$2"; shift 2 ;;
    --tile_size_km) TILE_SIZE_KM="$2"; shift 2 ;;
    --native_block_px) NATIVE_BLOCK_PX="$2"; shift 2 ;;
    --frequency) FREQUENCY="$2"; shift 2 ;;
    --pol_pair) POL_PAIR="$2"; shift 2 ;;
    --window) WINDOW="$2"; shift 2 ;;
    --noise_floor) NOISE_FLOOR="$2"; shift 2 ;;
    --min_coverage) MIN_COVERAGE="$2"; shift 2 ;;
    --use_exception_mask) USE_EXCEPTION_MASK="$2"; shift 2 ;;
    --lia) LIA="$2"; shift 2 ;;
    --lia_dem_scale) LIA_DEM_SCALE="$2"; shift 2 ;;
    --dem_dir) DEM_DIR="$2"; shift 2 ;;
    --nisar_dem_vrt) NISAR_DEM_VRT="$2"; shift 2 ;;
    --dem_cache_dir) DEM_CACHE_DIR="$2"; shift 2 ;;
    --earthdata_cookie) EARTHDATA_COOKIE="$2"; shift 2 ;;
    --gap_threshold) GAP_THRESHOLD="$2"; shift 2 ;;
    --bbox) BBOX="$2"; shift 2 ;;
    --keep_intermediate) KEEP_INTERMEDIATE="$2"; shift 2 ;;
    --script) SCRIPT_PATH="$2"; shift 2 ;;
    --python_exec) PYTHON_EXEC="$2"; shift 2 ;;
    --extra_args) EXTRA_ARGS="$2"; shift 2 ;;
    --grid_mode|--grid|--scale|--target_epsg|--grid_origin_x|--grid_origin_y)
      echo "$1 is not used here: the grid is always EASE2; choose --res 50 or 100" >&2; exit 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$INPUT_DIR" ]] || { echo "--input_dir is required" >&2; exit 1; }
[[ -n "$OUTPUT_DIR" ]] || { echo "--output_dir is required" >&2; exit 1; }
[[ -d "$INPUT_DIR" ]] || { echo "Input directory not found: $INPUT_DIR" >&2; exit 1; }
[[ -f "$SCRIPT_PATH" ]] || { echo "Python script not found: $SCRIPT_PATH" >&2; exit 1; }
command -v "$PYTHON_EXEC" >/dev/null 2>&1 || { echo "Python executable not found: $PYTHON_EXEC" >&2; exit 1; }
if [[ -n "$DEM_DIR" ]]; then
  [[ -d "$DEM_DIR" ]] || { echo "DEM directory not found: $DEM_DIR" >&2; exit 1; }
fi
[[ -n "$DEM_CACHE_DIR" ]] || DEM_CACHE_DIR="${INPUT_DIR}/DEM_Cache"

case "$RES" in 50|100) ;; *) echo "Invalid --res: $RES (use 50 or 100)" >&2; exit 1 ;; esac
case "$BACKSCATTER" in gamma0|sigma0) ;; *) echo "Invalid --backscatter: $BACKSCATTER" >&2; exit 1 ;; esac
[[ -z "$EASE2_MAT" || -f "$EASE2_MAT" ]] || { echo "EASE2 MAT not found: $EASE2_MAT" >&2; exit 1; }
case "$MODE" in nofilt|enhanced_lee|gamma_map|all) ;; *) echo "Invalid --mode: $MODE" >&2; exit 1 ;; esac
case "$MOSAIC" in yes|no) ;; *) echo "Invalid --mosaic: $MOSAIC" >&2; exit 1 ;; esac
case "$MOSAIC_ACROSS_MODES" in yes|no) ;; *) echo "Invalid --mosaic_across_modes: $MOSAIC_ACROSS_MODES" >&2; exit 1 ;; esac
case "$FREQUENCY" in frequencyA|frequencyB) ;; *) echo "Invalid --frequency: $FREQUENCY" >&2; exit 1 ;; esac
case "$POL_PAIR" in auto|HH_HV|VV_VH|RH_RV) ;; *) echo "Invalid --pol_pair: $POL_PAIR" >&2; exit 1 ;; esac
case "$KEEP_INTERMEDIATE" in yes|no) ;; *) echo "Invalid --keep_intermediate: $KEEP_INTERMEDIATE" >&2; exit 1 ;; esac
case "$USE_EXCEPTION_MASK" in yes|no) ;; *) echo "Invalid --use_exception_mask: $USE_EXCEPTION_MASK" >&2; exit 1 ;; esac
case "$LIA" in yes|no) ;; *) echo "Invalid --lia: $LIA" >&2; exit 1 ;; esac

mkdir -p "$OUTPUT_DIR"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"

launch_job() {
  local label="$1"
  local filter_name="$2"
  local destination="$3"
  mkdir -p "$destination"
  local log_file="${destination}/launcher_EASE2_${RES}m_${BACKSCATTER}_${label}_${TIMESTAMP}.log"

  local cmd=(
    "$PYTHON_EXEC" "$SCRIPT_PATH"
    --input_dir "$INPUT_DIR"
    --output_dir "$destination"
    --ease2_res "$RES"
    --backscatter "$BACKSCATTER"
    --frequency "$FREQUENCY"
    --pol_pair "$POL_PAIR"
    --filter "$filter_name"
    --filter_window "$WINDOW"
    --noise_floor "$NOISE_FLOOR"
    --min_coverage "$MIN_COVERAGE"
    --use_exception_mask "$USE_EXCEPTION_MASK"
    --lia "$LIA"
    --lia_dem_scale "$LIA_DEM_SCALE"
    --dem_cache_dir "$DEM_CACHE_DIR"
    --earthdata_cookie "$EARTHDATA_COOKIE"
    --mosaic "$MOSAIC"
    --mosaic_across_modes "$MOSAIC_ACROSS_MODES"
    --gap_threshold "$GAP_THRESHOLD"
    --workers "$WORKERS"
    --tile_workers "$TILE_WORKERS"
    --tile_size_km "$TILE_SIZE_KM"
    --native_block_px "$NATIVE_BLOCK_PX"
    --keep_intermediate "$KEEP_INTERMEDIATE"
  )

  [[ -z "$EASE2_MAT" ]] || cmd+=(--ease2_mat "$EASE2_MAT")

  if [[ -n "$DEM_DIR" ]]; then
    cmd+=(--dem_dir "$DEM_DIR")
  elif [[ -n "$NISAR_DEM_VRT" ]]; then
    cmd+=(--nisar_dem_vrt "$NISAR_DEM_VRT")
  fi

  if [[ -n "$BBOX" ]]; then
    read -r -a bbox_values <<< "$BBOX"
    [[ ${#bbox_values[@]} -eq 4 ]] || { echo "--bbox requires four values: W S E N" >&2; exit 1; }
    cmd+=(--bbox "${bbox_values[@]}")
  fi

  if [[ -n "$EXTRA_ARGS" ]]; then
    read -r -a extra_values <<< "$EXTRA_ARGS"
    cmd+=("${extra_values[@]}")
  fi

  nohup "${cmd[@]}" >"$log_file" 2>&1 &
  local pid=$!

  echo "Launched $label"
  echo "  Launcher version: $SCRIPT_VERSION"
  echo "  Processor       : $SCRIPT_PATH"
  echo "  PID             : $pid"
  echo "  Grid            : EASE-Grid 2.0, ${RES} m ($( [[ $RES == 50 ]] && echo M0p05 || echo M0p1 )), EPSG:6933"
  echo "  Backscatter     : $BACKSCATTER"
  echo "  Frequency       : $FREQUENCY"
  echo "  Filter          : $filter_name (window $WINDOW)"
  echo "  Mosaic          : $MOSAIC (across modes: $MOSAIC_ACROSS_MODES)"
  echo "  Exception mask  : $USE_EXCEPTION_MASK"
  echo "  LIA/incidence   : $LIA"
  if [[ "$LIA" == "yes" ]]; then
    echo "  Earthdata cookie: $EARTHDATA_COOKIE"
    echo "  DEM cache       : $DEM_CACHE_DIR"
  fi
  echo "  Output          : $destination"
  echo "  Log             : tail -f '$log_file'"
  echo
}

if [[ "$MODE" == "all" ]]; then
  echo "WARNING: --mode all starts three independent processors."
  launch_job nofilt none "${OUTPUT_DIR}/nofilt"
  launch_job enhanced_lee enhanced_lee "${OUTPUT_DIR}/enhanced_lee"
  launch_job gamma_map gamma_map "${OUTPUT_DIR}/gamma_map"
else
  case "$MODE" in
    nofilt) FILTER_NAME="none" ;;
    enhanced_lee) FILTER_NAME="enhanced_lee" ;;
    gamma_map) FILTER_NAME="gamma_map" ;;
  esac
  launch_job "$MODE" "$FILTER_NAME" "$OUTPUT_DIR"
fi
