# NISAR GCOV → EASE-Grid 2.0 (new standalone MAAP algorithm)

Independent repository for converting NISAR L2 GCOV HDF5 files to precisely
aligned **50 m (M0p05) or 100 m (M0p1)** EASE-Grid 2.0 (EPSG:6933) GeoTIFFs.

This package includes **unmodified copies** of the uploaded
`NISAR_GCOV_TO_EASE2.py` and `NISAR_GCOV_EASE2.sh` (2026-11-19). The Python
processor retains CRID-specific provisional and validated exception masks,
frequencyA/frequencyB handling, sigma0 conversion, filtering, LIA/incidence,
frame mosaicking, pixel-grid MAT companion files and no-overwrite/resume behavior.

## Contents

| File | Role |
|---|---|
| `NISAR_GCOV_TO_EASE2.py` | Original ~5,100-line processing algorithm, **unchanged** |
| `NISAR_GCOV_EASE2.sh` | Original asynchronous/`nohup` shell launcher for standalone HPC use, **unchanged** |
| `run_maap.sh` | **New** synchronous MAAP DPS run command, with file staging and `output/` collection |
| `build-env.sh` | **New** MAAP build command (conda environment) |
| `environment.yml` | Dependencies: Python, HDF5, NumPy, SciPy, Rasterio, GDAL and Matplotlib |
| `algorithm_config.yaml` | Reference MAAP algorithm metadata and named inputs |
| `MAAP_REGISTRATION.md` | Complete registration form values, limitations and GitHub instructions |
| `tests/test_maap_wrapper.sh` | Wrapper smoke test using a stub processor (no NISAR data required) |

## Defaults

- Grid: EASE2 **100 m** (actual lattice pixel size 100.089502334956 m)
- Radiometry: **gamma0**, stored in **linear power**; quicklooks use dB
- Frequency: `frequencyA`; polarization pair: `auto`
- Filter: `nofilt`; optional enhanced Lee and Gamma MAP filters
- Mosaic: `yes` (compatible frames only; one-file DPS jobs do not mosaic separate jobs)
- Product maturity: provisional vs validated `inputDataExceptionMask` rules
- Local incidence angle: `no` (option available after DEM credentials/setup)
- No geographic clipping unless `--bbox` is explicitly provided

## MAAP DPS

Use **`run_maap.sh`**, not the original background launcher:

```bash
bash run_maap.sh --input_file input/NISAR_L2_PR_GCOV_....h5 \
  --ease2_res 100 --backscatter gamma0 --mode nofilt
```

The runner auto-creates `output/`, executes the processor in the foreground,
returns errors to MAAP, and writes logs and GeoTIFFs into the directory that
DPS collects. By default it runs with the conda environment built by
`build-env.sh`.

For local use of the original cluster launcher:

```bash
bash NISAR_GCOV_EASE2.sh \
  --input_dir /path/to/gcov --output_dir /path/to/results \
  --res 50 --backscatter sigma0 --mode enhanced_lee
```

Read [MAAP_REGISTRATION.md](MAAP_REGISTRATION.md) for the exact **new**
GitHub repo name and registration fields, the one-granule-per-job default,
and `--lia yes` DEM-access requirements.

## Checks performed on this package

```bash
bash -n run_maap.sh build-env.sh NISAR_GCOV_EASE2.sh
python -m py_compile NISAR_GCOV_TO_EASE2.py
bash tests/test_maap_wrapper.sh
```

These checks validate syntax and the MAAP wrapper's input resolution / output
layout. They **do not validate numerical GCOV processing** without real HDF5
granules; a MAAP staging and full processing test is still required.
