# Register this NEW algorithm on NASA MAAP (DPS)

**Important:** Register with a **new, unique algorithm name**. Reusing an existing
MAAP algorithm name can replace the previous algorithm, even if the repository
is different. The original `NISAR_GCOV_TO_TIFF` repo must not be changed.

After creating a new **public** GitHub repository with these files at its root:

| Registration field | Value |
|---|---|
| Repository URL | `https://github.com/preet0lal/nisar-gcov-ease2-maap.git` |
| Repository Branch | `main` |
| Run Command | `nisar-gcov-ease2-maap/run_maap.sh` |
| Build Command | `nisar-gcov-ease2-maap/build-env.sh` |
| Algorithm Name | `nisar-gcov-ease2-v1` |
| Algorithm Description | `NISAR L2 GCOV to EPSG:6933 EASE2 50/100 m linear gamma0/sigma0 GeoTIFF, CRID-aware masks, filtering, mosaics and optional LIA.` |
| Disk Space (GB) | `200` (starting point; adapt to frame size, resolution and filtering) |
| Resource Allocation | `maap-dps-worker-32vcpu-64gb` **if present in your dropdown**; otherwise select an appropriate available queue |
| Container URL | `mas.maap-project.org/root/maap-workspaces/custom_images/maap_base:v6.0.0` |

MAAP may display the code path relative to a clone directory rather than the
repository root. The example above follows MAAP's current documented
`repository-name/path/to/script.sh` convention. If the build fails with a file
not found, check the clone path printed in the registration log.

## Inputs to add in the MAAP form

### File Inputs

Create one **required** File Input with name `input_file`. Its value at job
submission should be a MAAP-accessible HDF5 file URI (HTTPS or S3), not a path
on your personal workstation or an inaccessible MAAP notebook filesystem.
DPS stages the file before running the script. The wrapper handles the staged
local file or resolves a URI to the staged basename.

### Optional named parameters

Use the optional named arguments in the MAAP interface if that interface
supports them. Some registration screens call these *Positional Inputs*, but
the actual run script expects `--name value`, not unnamed positional values.

| Parameter | Default | Allowed values / purpose |
|---|---|---|
| `ease2_res` | `100` | `50`, `100` |
| `backscatter` | `gamma0` | `gamma0`, `sigma0` (linear power) |
| `frequency` | `frequencyA` | `frequencyA`, `frequencyB` |
| `pol_pair` | `auto` | `auto`, `HH_HV`, `VV_VH`, `RH_RV` |
| `mode` | `nofilt` | `nofilt`, `enhanced_lee`, `gamma_map`, `all` |
| `filter_window` | `7` | odd integer ≥3 |
| `mosaic` | `yes` | `yes`, `no` |
| `use_exception_mask` | `yes` | `yes`, `no`; preserves product-maturity logic |
| `lia` | `no` | `yes`, `no`; see credential warning below |
| `workers` | `1` | number of acquisition-group processes |
| `tile_workers` | `4` | threads per processing group |

**Simplest registration:** add only the `input_file` File Input and no other
inputs. The wrapper then uses 100 m, gamma0, frequencyA, no filter, and CRID-aware
quality masks. To change settings through the MAAP Submit Job form, explicitly
register the matching named parameters; otherwise run the wrapper from a
terminal using additional flags.

## Data layout and limitations

- The MAAP runner writes to the working directory's `output/` in the **foreground**.
  Output TIFF, PNG, MAT, and logs are uploaded by DPS upon successful completion.
- The untouched original `NISAR_GCOV_EASE2.sh` is for cluster/terminal use and
  starts a detached `nohup` job. Do **not** use it as the MAAP Run Command.
- One File Input means **one GCOV granule per MAAP job**. Compatible adjacent
  frames can be mosaicked only if multiple granules are staged in the **same**
  job and supplied to `--input_dir`. A notebook's local file path cannot be
  assumed to exist on a remote DPS worker.
- The processor preserves the existing no-overwrite semantics. `--bbox` is
  **not supplied** by default, so SAR is not clipped to an AOI.
- For `--lia yes`, the processor needs an accessible NISAR DEM: either HTTPS
  with valid Earthdata credentials within the worker or a staged local DEM.
  Neither `~/.netrc` nor a local DEM cache from your notebook is transferred
  automatically to DPS. Start with `lia=no` to validate the registration.
- EASE2 50 m can require substantial temporary storage. Choose sufficient disk
  and RAM based on GCOV frame size; `200 GB` is a suggested starting point,
  not a benchmarked requirement.

## Test in a MAAP terminal before registration

```bash
bash build-env.sh
bash tests/test_maap_wrapper.sh
bash run_maap.sh --help
# For a REAL local test, with existing GCOV files:
mkdir -p output
bash run_maap.sh --input_dir /path/to/gcov_frames \
  --ease2_res 100 --backscatter gamma0 --mode nofilt
```

## Push to a separate GitHub repository

First, create a **new empty, public** repo named `nisar-gcov-ease2-maap`
at https://github.com/new under `preet0lal` (do not select any existing repo).
When using this folder from the ZIP:

```bash
cd nisar-gcov-ease2-maap
git init -b main
git add .
git commit -m "Add standalone NISAR GCOV EASE2 MAAP DPS algorithm"
git remote add origin https://github.com/preet0lal/nisar-gcov-ease2-maap.git
git push -u origin main
```

Do **not** add this as a new branch inside the old NISAR GCOV repository.
