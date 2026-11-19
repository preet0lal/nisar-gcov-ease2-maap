# NISAR GCOV → EASE-Grid 2.0 (50 m / 100 m)

Convert **NASA–ISRO SAR (NISAR) Level-2 Geocoded Covariance (GCOV)** HDF5 granules into georeferenced **EASE-Grid 2.0 GeoTIFFs** at nominal **50 m or 100 m**, with selectable **gamma-naught (γ⁰)** or **sigma-naught (σ⁰)** backscatter. Designed for processing on a Linux workstation/HPC cluster or as a **NASA MAAP Data Processing System (DPS)** algorithm.

> **Separate project.** This is a NEW repository. It does not change the older GCOV-to-TIFF repository. `NISAR_GCOV_TO_EASE2.py` and `NISAR_GCOV_EASE2.sh` are the original, unchanged EASE2 scripts supplied on 2026-11-19; the MAAP integration is added alongside them.

## At a glance

| Item | Details |
|---|---|
| Input | NISAR L2 GCOV `.h5` / `.hdf5` files, using standard NISAR filenames |
| Output projection | **EPSG:6933**, global EASE-Grid 2.0 cylindrical equal-area |
| Output resolution | Nominal **100 m (M0p1)**, or **50 m (M0p05)** |
| Backscatter | `gamma0` (default) or `sigma0`, **linear power**, not dB |
| Polarizations | Auto-select supported co/cross pair; or choose `HH_HV`, `VV_VH`, `RH_RV` |
| Frequency | `frequencyA` (default) or `frequencyB` |
| Speckle | No filter (default), enhanced Lee, or Gamma MAP |
| Quality screening | Subswath, looks, finite covariance, optional CRID-aware exception mask (default on) |
| Optional outputs | Incidence angle, terrain local incidence angle (LIA), quicklook, EASE2 row/column `.mat` companion |
| Frames | Compatible adjacent frames can be mosaicked; different dates/tracks/passes are kept separate |
| Resume | Existing complete output sets are skipped; conflicting/partial sets are **not overwritten** |

## Which script should I run?

| File | Use |
|---|---|
| `NISAR_GCOV_TO_EASE2.py` | The underlying scientific processor; accepts CLI arguments and runs in the foreground |
| `NISAR_GCOV_EASE2.sh` | Original **local/HPC** launcher; starts background `nohup` jobs and prints logs/PIDs |
| **`run_maap.sh`** | **MAAP DPS entry point**; runs in the **foreground** and writes to DPS `output/` |
| `build-env.sh` / `environment.yml` | Creates the MAAP-compatible `nisar_gcov_ease2` Conda environment and checks base-image S3 helpers |
| `maap_s3.py` | Obtains temporary Earthdata S3 credentials on MAAP, and downloads an `s3://` GCOV input if DPS did not stage it |
| `MAAP_REGISTRATION.md` | Registration form settings and operational caveats |
| `algorithm_config.yaml` | Human-readable example registration configuration; **not** itself automatically imported into MAAP |
| `tests/test_maap_wrapper.sh` | Fast wrapper smoke test without NISAR data |

## 1. Install and run locally

Prerequisites: Linux, a Conda installation with `conda` in `PATH`, and disk space for temporary GeoTIFFs. GDAL utilities are needed when generating LIA from a DEM.

```bash
git clone https://github.com/preet0lal/nisar-gcov-ease2-maap.git
cd nisar-gcov-ease2-maap
bash build-env.sh
conda activate nisar_gcov_ease2
python NISAR_GCOV_TO_EASE2.py --help
```

### Example A: 100 m gamma0 (default, no filter)

```bash
python NISAR_GCOV_TO_EASE2.py \
  --input_dir /data/NISAR/GCOV \
  --output_dir /data/NISAR/EASE2_100m_gamma0 \
  --ease2_res 100 --backscatter gamma0 \
  --filter none --mosaic yes --use_exception_mask yes
```

### Example B: 50 m sigma0 with enhanced Lee filtering

```bash
python NISAR_GCOV_TO_EASE2.py \
  --input_dir /data/NISAR/GCOV \
  --output_dir /data/NISAR/EASE2_50m_sigma0 \
  --ease2_res 50 --backscatter sigma0 \
  --filter enhanced_lee --filter_window 7
```

### Example C: Original background launcher on HPC

```bash
bash NISAR_GCOV_EASE2.sh \
  --input_dir /data/NISAR/GCOV \
  --output_dir /data/NISAR/EASE2_100m \
  --res 100 --backscatter gamma0 --mode nofilt \
  --workers 2 --tile_workers 4
```

The **original** launcher immediately returns and prints a PID and a `tail -f` command. It is **not suitable as the MAAP Run Command** because DPS must monitor the foreground process to know when it has finished.

## 2. Grid definition and geospatial alignment

Both output grids use the exact **global EASE2 lattice**, not a new grid origin calculated for each scene. Processing crops this lattice to the coverage of the selected GCOV frame(s) while retaining globally consistent cell boundaries.

| Property | 100 m | 50 m |
|---|---:|---:|
| EASE2 grid ID | M0p1 | M0p05 |
| Exact cell spacing (m) | 100.089502334956 | 50.044751167478 |
| Full global columns | 347,040 | 694,080 |
| Full global rows | 146,160 | 292,320 |
| Cells per 1 km EASE2 M01 cell | 10 × 10 | 20 × 20 |

One 100 m cell is exactly **2 × 2** 50 m cells. The pixel-center row/column coordinates and their parent 1 km M01 identifiers are saved in a companion `*_EASE2_grid.mat` file. The MATLAB `.mat` stores **1-D row and column vectors**, not a huge global 2-D array. Row/column IDs are **zero-based** by default, with one-based versions included for MATLAB.

The optional `--ease2_mat` checks a reference latitude/longitude MAT against the selected grid; it is **not needed** for standard processing.

**Full-frame behavior:** By default, the program does **not** clip to any AOI. Each selected acquisition contributes its full frame extent. Supplying `--bbox 'W S E N'` explicitly **does clip** the output to that geographic box. To *select frames that intersect an AOI yet keep full frames*, select those GCOV files upstream and **omit** `--bbox` during conversion.

## 3. Radiometry and quality masking

### Gamma0 versus sigma0

The GCOV product's stored normalization is read from its metadata. If conversion is necessary, it is applied **at native resolution, before speckle filtering and EASE2 aggregation**:

```text
sigma0_linear = gamma0_linear * rtcGammaToSigmaFactor
gamma0_linear = sigma0_linear / rtcGammaToSigmaFactor
```

If the file already uses the requested convention, it is left unchanged. If the metadata field is missing, the processor warns and assumes the source is gamma0. **The exported SAR GeoTIFF values are linear power.** The PNG quicklook converts positive values to dB for visual inspection only: `dB = 10 log10(linear)`.

### Provisional versus validated products

With `--use_exception_mask yes` (default), the processor applies different `inputDataExceptionMask` rules depending on the **Composite Release ID (CRID)** and the mask's stored type:

| Product | Mask format | Retain pixel when… |
|---|---|---|
| **Provisional**, CRID number below 5030 | `uint8` | Exception-mask value is **0** |
| **Validated**, CRID number 5030 or above | `uint16` | Not fill (typically 65535), anomaly low-byte bits are zero, **and both selected polarization validity bits are set** |

For validated HH/HV, the required bits are HH=`256` and HV=`512`, so `768` is a typical good value—not zero. For older provisional **frequencyB** files whose exception layer has frequencyA dimensions, the layer is reduced by a *k×k* all-valid rule (typically **8×8**): if any fine pixel is flagged, its coarse frequencyB pixel is rejected.

The processor also rejects pixels failing the GCOV subswath mask (`1..numberOfSubSwaths`), `numberOfLooks > 0`, or finite/nonnegative covariance checks. `--min_coverage 0.90` requires at least 90% valid source-area coverage per output pixel. With `--use_exception_mask no`, exception-mask exclusions are disabled **but all other validity requirements remain**.

## 4. Processing, filtering, and mosaics

1. Find GCOV files and read their track, pass, frame, acquisition date, mode, CRID, frequency, and polarization metadata.
2. Apply native quality masks and any gamma0/sigma0 conversion; optionally filter native-grid backscatter (`none`, `enhanced_lee`, `gamma_map`).
3. Reproject/area-average **linear** values onto the globally aligned EASE2 grid; enforce minimum valid-area coverage.
4. Optionally combine adjacent, compatible frames with a per-pixel mean and save a frame-count image. Different dates, tracks, or passes are not mixed. Different bandwidth modes are kept separate by default (`--mosaic_across_modes no`).
5. Write output GeoTIFFs, a quicklook, and the EASE2 MATLAB grid companion; optionally write incidence and LIA GeoTIFFs.

**Important:** `--mosaic yes` only combines compatible frames **present in the same processing invocation**. If MAAP submits one GCOV granule per DPS job, separate jobs will **not** automatically mosaic across one another. To mosaic adjacent frames, stage the frames together and run the processor with one input directory.

**Resolution caution:** `frequencyB` can be coarser than 50 m. Exporting it at nominal 50 m produces an aligned **oversampled** grid; it does **not** create additional spatial information.

## 5. Outputs

Files are grouped by cycle, track, pass, compatible frame list, mode, polarization, frequency, date, CRID, EASE2 resolution, radiometry, and filter. The filename stem follows this pattern:

```text
NISAR_L_GCOV_c<cycle>_t<track>_<A|D>_f<frames>_<mode>_<pol>_<frequency>_<date>_<CRID>_EASE2_<50|100>m_<gamma0|sigma0>_<filter>_
```

| File suffix | Contents |
|---|---|
| `_*co-pol-term*.tif` (e.g. `_HHHH.tif`) | Co-polarized backscatter, linear power, Float32 |
| `_*cross-pol-term*.tif` (e.g. `_HVHV.tif`) | Cross-polarized backscatter, linear power, Float32 |
| `_frame_count.tif` | Count of valid contributing frames per pixel |
| `_quicklook.png` | SAR visualization in dB; geometry panels when available |
| `_EASE2_grid.mat` | Global row/column IDs, parent 1 km IDs, x/y, lon/lat pixel-center vectors |
| `_incidence_angle_deg.tif` | Optional ellipsoid-referenced incidence angle in degrees |
| `_local_incidence_angle_deg.tif` | Optional DEM-derived local incidence angle in degrees |
| `processing.log` | Processing diagnostics and summary |

GeoTIFF outputs have EPSG:6933 georeferencing and compressed, tiled storage. Default SAR NoData is `-9999`. Temporary native-frame files are removed by default unless `--keep_intermediate yes` is selected.

## 6. Register as a NEW MAAP DPS algorithm

Open **MAAP → Register Algorithm**, with the code committed in this public repository. Do **not** select or overwrite the previous GCOV algorithm. Use the new algorithm name shown below.

| MAAP field | Value |
|---|---|
| Repository URL | `https://github.com/preet0lal/nisar-gcov-ease2-maap.git` |
| Repository Branch | `main` |
| Run Command | `nisar-gcov-ease2-maap/run_maap.sh` |
| Build Command | `nisar-gcov-ease2-maap/build-env.sh` |
| Algorithm Name | `nisar-gcov-ease2-v1` |
| Disk Space | **200 GB** (starting estimate; size for actual granules) |
| Resource Allocation | Choose an available worker with enough CPU/RAM; `maap-dps-worker-32vcpu-64gb` if offered |
| Container URL | `mas.maap-project.org/root/maap-workspaces/custom_images/maap_base:v6.0.0` |

**Required File Input:** Add one named `input_file` (a MAAP-accessible GCOV HDF5 URI). At job submission, DPS can stage it to the worker; `run_maap.sh` accepts the staged path, resolves the original URL by basename in known input directories, or attempts a direct `s3://` download using temporary Earthdata S3 credentials via `maap_s3.py`. Start with one granule and default parameters. The runner writes to the job's `output/` and executes **synchronously**, so MAAP can collect completed results and report failures.

### MAAP positional input order

The **current GitHub `run_maap.sh` supports both** named flags (`--ease2_res 50`) **and bare positional values** supplied by MAAP DPS. When adding *Positional Inputs* in the registration UI, use **exactly this order after the File Input**; `none`, `-`, `null`, or an empty value means “keep the wrapper default”:

| Position | MAAP input name | Default | Allowed values |
|---:|---|---|---|
| 1 | `input_file` (**File Input**) | Required | One GCOV HDF5 URI, e.g. `s3://…/NISAR_…h5` |
| 2 | `ease2_res` | `100` | `50`, `100` |
| 3 | `mode` | `nofilt` | `nofilt`, `enhanced_lee`, `gamma_map`, `all` |
| 4 | `filter_window` | `7` | Odd integer ≥ 3 |
| 5 | `lia` | `no` | `yes`, `no` |
| 6 | `frequency` | `frequencyA` | `frequencyA`, `frequencyB` |
| 7 | `backscatter` | `gamma0` | `gamma0`, `sigma0` |
| 8 | `clip_to_aoi` | `no` when no bbox is supplied | `yes` to crop; `no` for full frames |
| 9–12 | `bbox_w`, `bbox_s`, `bbox_e`, `bbox_n` | Unset | West, south, east, north, decimal degrees; only needed when clipping |

Example **bare positional** call (50 m, enhanced Lee, sigma0, full frame):

```bash
bash run_maap.sh s3://bucket/path/NISAR_L2_PR_GCOV_EXAMPLE.h5 \
  50 enhanced_lee 7 no frequencyA sigma0 no
```

This sample URI is illustrative, not an existing granule. For a first registration, you can add **only the required File Input** and use defaults. The named-flag interface remains available for manual runs.

Full registration details: [MAAP_REGISTRATION.md](MAAP_REGISTRATION.md).

### Test the MAAP wrapper in a terminal

```bash
bash tests/test_maap_wrapper.sh
# Direct MAAP-runner example with a locally accessible file:
bash run_maap.sh --input_file /data/NISAR/GCOV/example.h5 \
  --ease2_res 100 --backscatter gamma0 --mode nofilt
```

The example name `example.h5` is illustrative; use an actual **standard-named** NISAR L2 GCOV granule. Running `run_maap.sh` from a MAAP job directory creates that directory's `output/`.

## 7. Generate incidence and local incidence angles

LIA is **disabled by default** because it requires terrain data and additional processing. To enable it locally with existing DEM tiles:

```bash
python NISAR_GCOV_TO_EASE2.py \
  --input_dir /data/NISAR/GCOV --output_dir /data/NISAR/EASE2_LIA \
  --ease2_res 100 --lia yes --lia_dem_scale 30 \
  --dem_dir /data/NISAR/DEM_tiles
```

Alternatively, the script can read the official NISAR DEM VRT over HTTPS. On MAAP, `run_maap.sh` first attempts to obtain **temporary Earthdata S3 credentials** through `maap_s3.py` and use the NISAR DEM VRT via `/vsis3`; if this fails it falls back to the processor's default DEM URL, which may require authenticated HTTPS access. The environment also needs GDAL tools `gdalwarp` / `gdalbuildvrt`. MAAP DPS workers **do not automatically inherit** credentials or DEM caches from your laptop/notebook. For a first MAAP test, leave `--lia no`.

The resulting angle rasters use the **same EASE2 grid and SAR-valid coverage** as the backscatter outputs. They are angles in **degrees**, not backscatter quantities.

## 8. Common configuration options

| Flag (Python / `run_maap.sh`) | Default | Notes |
|---|---|---|
| `--input_dir` / `--input_file` | Required | Python processor accepts a **directory**; MAAP wrapper can accept one staged file |
| `--ease2_res` | `100` | `50` or `100` |
| `--backscatter` | `gamma0` | `gamma0` / `sigma0` |
| `--frequency` | `frequencyA` | `frequencyA` / `frequencyB` |
| `--pol_pair` | `auto` | `auto`, `HH_HV`, `VV_VH`, `RH_RV` |
| `--filter` (Python) / `--mode` (wrappers) | `none` / `nofilt` | Enhanced Lee or Gamma MAP; wrapper `all` runs each filter separately |
| `--filter_window` | `7` | Odd integer ≥ 3 |
| `--noise_floor` | `0` | Linear-power threshold; zero disables |
| `--min_coverage` | `0.90` | Native valid-area fraction per output cell |
| `--use_exception_mask` | `yes` | CRID-aware quality rule |
| `--mosaic` | `yes` | Compatible adjacent input frames |
| `--mosaic_across_modes` | `no` | Keep bandwidth modes separate unless explicitly enabled |
| `--gap_threshold` | `50` | Maximum inter-frame gap, km |
| `--lia` | `no` | Optional incidence + local incidence angle |
| `--lia_dem_scale` | `30` | Internal DEM processing spacing, m |
| `--workers` | `1` | Parallel acquisition groups |
| `--tile_workers` | `4` | Warp threads per group |
| `--tile_size_km` | `50` | Output window size for processing |
| `--native_block_px` | `1024` | Native HDF5 processing block dimension |
| `--bbox` | Unset | Geographic `W S E N`; wrapper **clips** if bbox supplied and `clip_to_aoi` is not `no` |
| `--clip_to_aoi` (MAAP wrapper) | `no` without bbox | `yes` uses bbox; `no` keeps complete frames even when bbox values were supplied |
| `--keep_intermediate` | `no` | Retain temporary native products |
| `--ease2_mat` | Unset | Optional grid coordinate validation |
| `--dem_dir`, `--nisar_dem_vrt`, `--dem_cache_dir` | Optional | Terrain source/cache for LIA |

Run `python NISAR_GCOV_TO_EASE2.py --help` and `bash run_maap.sh --help` for the supported interfaces. The underlying Python processor itself has `--bbox` but does **not** have a `--clip_to_aoi` option; that switch belongs to the MAAP wrapper. `NISAR_GCOV_EASE2.sh` uses `--res` and `--window` instead of the Python flag names.

## 9. Verification and known limitations

```bash
bash -n NISAR_GCOV_EASE2.sh run_maap.sh build-env.sh
python -m py_compile NISAR_GCOV_TO_EASE2.py
bash tests/test_maap_wrapper.sh
```

The smoke test uses a **dummy processor**, so it verifies MAAP file staging, arguments, foreground execution, and output-directory routing—not the scientific results. A real GCOV HDF5 smoke test on the MAAP DPS worker is still necessary before relying on production outputs.

- **No real-data numerical validation is claimed here.** Validate projection, masks, radiometry, and output values on a known granule.
- **No implicit scene selection by AOI.** Select the input granules yourself. `--bbox` is strictly an output clip.
- **No mosaicking across independent DPS jobs.** Place compatible frames in the same invocation to mosaic them.
- **No data beyond source spatial information.** A finer output grid does not make coarse native posting more informative.
- **No-overwrite by design.** To rerun a changed configuration that reuses the same output stem, use a clean output directory rather than expecting old files to be replaced.
- **Disk and RAM can be substantial**, especially for 50 m data, filtering, multiple frames, or LIA. MAAP queue/disk suggestions are starting points, not measured minimum requirements.

## 10. Troubleshooting

| Symptom | What to check |
|---|---|
| MAAP says input missing | Check DPS staging to `input/`, `inputs/`, or `/inputs`; for an un-staged `s3://` URI verify that MAAP can obtain temporary Earthdata S3 credentials with `maap_s3.py` |
| MAAP exits while background task continues | Registered **`run_maap.sh`**, not the `nohup` HPC launcher |
| Missing `conda` or Python modules | Confirm `build-env.sh` ran with a container providing Conda; inspect build logs |
| All pixels masked | Check CRID, `inputDataExceptionMask` dtype/bits, native subswath mask, and `--min_coverage`; validated good pixels need not have mask value zero |
| `--lia yes` fails | Check authenticated DEM access, GDAL utilities, local DEM coverage, and cache writability |
| Existing acquisition skipped | Check processing log for preexisting partial files or metadata mismatch; existing outputs are deliberately protected |
| Need complete intersecting frames | **Omit `--bbox`** and provide the GCOV files intersecting your AOI |

## Provenance

EASE2-only processing scripts provided **2026-11-19**; MAAP DPS wrapper and documentation added in this independent repository. The source-code checksums are available in [`SOURCE_SHA256.txt`](SOURCE_SHA256.txt). There is no stated license for these source files; contact the repository owner before assuming permission to redistribute or reuse outside applicable terms.
