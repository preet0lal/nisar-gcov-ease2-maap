# Register NISAR GCOV EASE2 on MAAP DPS

This guide belongs to the **new and independent** repository [`preet0lal/nisar-gcov-ease2-maap`](https://github.com/preet0lal/nisar-gcov-ease2-maap). It does **not** replace or modify the older `NISAR_GCOV_TO_TIFF` repository or algorithm. Do not reuse the existing MAAP algorithm name.

## 1. MAAP registration form

| MAAP field | Enter this |
|---|---|
| Repository URL | `https://github.com/preet0lal/nisar-gcov-ease2-maap.git` |
| Repository Branch | `main` |
| Run Command | `nisar-gcov-ease2-maap/run_maap.sh` |
| Build Command | `nisar-gcov-ease2-maap/build-env.sh` |
| Algorithm Name | `nisar-gcov-ease2-v1` |
| Algorithm Description | `Convert NISAR L2 GCOV HDF5 to EASE-Grid 2.0 50/100 m linear gamma0/sigma0 GeoTIFFs with CRID-aware masks, filtering, mosaics and optional incidence/LIA.` |
| Disk Space (GB) | `200` to start; tune for actual scene size and options |
| Resource Allocation | Choose a queue with sufficient CPU/RAM; `maap-dps-worker-32vcpu-64gb` if available |
| Container URL | `mas.maap-project.org/root/maap-workspaces/custom_images/maap_base:v6.0.0` |

The run/build paths above follow MAAP's `repository-name/script.sh` convention. If the MAAP build log shows a different checkout path, use the path shown there.

## 2. File Input and Positional Inputs

Create a required **File Input** named `input_file`. At submission, use a MAAP-accessible URI to **one standard-named NISAR L2 GCOV HDF5 granule** (usually `.h5`). MAAP may stage the file automatically. If an `s3://` input is not staged, `run_maap.sh` attempts to download it using temporary Earthdata S3 credentials from MAAP through `maap_s3.py`.

**Important: `run_maap.sh` accepts raw positional inputs**, in the following exact order. Register optional *Positional Inputs* in this order **after** `input_file` if you need to control them from the MAAP submission form:

| Overall position | Input name | Default | What to enter |
|---:|---|---|---|
| 1 | `input_file` (**File Input**) | Required | One MAAP-accessible HDF5 URI |
| 2 | `ease2_res` | `100` | `50` or `100` |
| 3 | `mode` | `nofilt` | `nofilt`, `enhanced_lee`, `gamma_map`, `all` |
| 4 | `filter_window` | `7` | Odd integer ≥ 3 |
| 5 | `lia` | `no` | `yes` or `no` |
| 6 | `frequency` | `frequencyA` | `frequencyA` or `frequencyB` |
| 7 | `backscatter` | `gamma0` | `gamma0` or `sigma0`, linear power |
| 8 | `clip_to_aoi` | `no` if bbox not supplied | `no` to preserve entire frame; `yes` to crop |
| 9 | `bbox_w` | Unset | West longitude, decimal degrees |
| 10 | `bbox_s` | Unset | South latitude, decimal degrees |
| 11 | `bbox_e` | Unset | East longitude, decimal degrees |
| 12 | `bbox_n` | Unset | North latitude, decimal degrees |

In the **positional** interface, `none`, `-`, `null`, or empty values retain defaults. To crop, provide all four bounds in EPSG:4326, with west < east and south < north. To **keep whole GCOV frames even when they intersect only part of an AOI**, set `clip_to_aoi=no` (or omit bbox arguments entirely).

### Example 1 — defaults

Register only `input_file`. The wrapper will use **100 m**, **gamma0**, **nofilt**, **frequencyA**, **LIA off**, default quality masking, and **full-frame output**.

### Example 2 — 50 m sigma0, enhanced Lee, full frame

The equivalent manual positional command is:

```bash
bash run_maap.sh s3://bucket/path/standard_named_NISAR_GCOV.h5 \
  50 enhanced_lee 7 no frequencyA sigma0 no
```

The URI above is an illustrative placeholder; use the real granule URL. The order of registered optional inputs must match the order in this command.

### Example 3 — enable an explicit output clip

```bash
bash run_maap.sh s3://bucket/path/standard_named_NISAR_GCOV.h5 \
  100 nofilt 7 no frequencyA gamma0 yes 95 -11 141 6
```

This example crops output to west=95°, south=-11°, east=141°, north=6°; it is **not** a requirement of the algorithm.

## 3. Named flag interface for debugging and local runs

The same runner accepts named flags. Unlike the original HPC launcher, it runs in the foreground and writes to the MAAP job working directory's `output/`:

```bash
bash run_maap.sh \
  --input_file /path/to/actual_NISAR_GCOV.h5 \
  --ease2_res 100 --mode nofilt --backscatter gamma0 \
  --lia no --clip_to_aoi no
```

The Python processor itself uses `--input_dir` and `--output_dir`, and `--filter none` instead of wrapper `--mode nofilt`; see the main README.

## 4. Important processing and infrastructure notes

- **MAAP must run `run_maap.sh`**, not `NISAR_GCOV_EASE2.sh`. The latter launches detached `nohup` processes and is for standalone HPC use.
- One `input_file` means one granule per DPS job. `--mosaic yes` only mosaics compatible frames **within one invocation**, not across separate MAAP jobs.
- If a granule is supplied by S3 and not staged, the runner requests temporary Earthdata S3 credentials via `maap_s3.py` and downloads the file. Access depends on MAAP authentication, bucket permissions, and availability of MAAP's base Python helper dependencies (`maap-py`, `s3fs`).
- For `lia=yes`, the runner attempts to use temporary Earthdata S3 credentials and the NISAR DEM via `/vsis3`, otherwise falls back to its default DEM URL. GDAL access and DEM availability are required. A local laptop's `~/.netrc`, cookies, or cache do not automatically follow a DPS job. Start with `lia=no` for first deployment testing.
- Outputs are saved under `output/`: linear backscatter GeoTIFFs, a frame-count GeoTIFF, quicklook PNG, EASE2 grid MAT, logs, and optional geometry rasters.
- **No overwrite:** complete existing acquisitions are skipped; incompatible or partial results are not silently replaced.
- 200 GB and a 32-vCPU/64-GB queue are **initial suggestions**, not validated minimum requirements. 50 m processing and LIA can require much more scratch space.

## 5. Test before production

```bash
bash -n run_maap.sh build-env.sh NISAR_GCOV_EASE2.sh
python -m py_compile NISAR_GCOV_TO_EASE2.py
bash tests/test_maap_wrapper.sh
bash run_maap.sh --help
```

The provided wrapper smoke test uses a dummy processor. It checks argument handling and output routing but is **not a real-data numerical validation**. Run one known GCOV acquisition end to end on MAAP, check the `processing.log`, confirm EPSG:6933 and the expected masks/normalization, then expand to batches.


## NEW: one DPS job for every ISMN frame (manifest batch mode)

The updated `run_maap.sh` now accepts either an individual `.h5`/`.hdf5` file (old behavior) or **one `.csv` manifest file** (new behavior). For CSV, the job sequentially processes every HDF5 listed in the manifest without creating any child DPS jobs. No scientific processing code was changed.

**First re-register/rebuild** `nisar-gcov-ease2-v1` from the latest `main` commit. Pushing GitHub code alone does not update an already-registered MAAP container. Keep the original registration's ordered inputs: `input_file, ease2_res, mode, filter_window, lia, frequency, backscatter, clip_to_aoi`. Re-registering this algorithm name replaces its previous MAAP build. This leaves the older, separate NISAR-to-TIFF GitHub repository untouched.

Open [`notebooks/ISMN_NISAR_GCOV_EASE2_50m_sigma0_ONE_JOB.ipynb`](notebooks/ISMN_NISAR_GCOV_EASE2_50m_sigma0_ONE_JOB.ipynb) in MAAP JupyterLab alongside `ismn_gt50_coords_points_june.csv`. Run discovery, inspect the manifest, then enable its confirmation flags. It uploads the manifest into your authorized workspace S3 prefix, and makes **one** `maap.submitJob` call. It uses 50 m EPSG:6933, linear sigma0, enhanced Lee 7x7, frequencyA, no AOI clipping and no LIA. The worker downloads each HDF5 using NISAR Earthdata credentials, processes full-frame outputs, deletes the HDF5 and proceeds to the next.

Batch output files include `single_dps_job_scene_ledger.csv`, `single_dps_job_summary.json`, the manifest copy and individual SAR GeoTIFF/PNG/MAT outputs under the **same** job's `output/` folder. A nonzero result indicates at least one failed scene. Batch mode saves separate scene outputs and does **not** mosaic adjacent frames across the manifest because it processes one frame at a time.

**Limits:** a typical DPS job has a **24-hour maximum** and a finite disk allocation; all 50 m scenes may not finish in a single attempt. Submit only after reviewing the number of matched scenes and the estimated processing time. DPS generally publishes output upon job completion; a failed or timed-out job may have incomplete or triaged outputs.
