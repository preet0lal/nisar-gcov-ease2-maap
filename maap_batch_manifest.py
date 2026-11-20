#!/usr/bin/env python3
"""Process a GCOV CSV manifest in ONE MAAP DPS job (sequential, full-frame).

Designed for the existing nisar-gcov-ease2-v1 run_maap.sh wrapper. The MAAP
File Input points to a CSV manifest, not an individual HDF5. Each listed NISAR
GCOV S3 URI is downloaded in turn using MAAP Earthdata credentials, processed
by run_maap.sh in the foreground, then deleted to limit scratch usage.

IMPORTANT: a single DPS job can hit runtime, disk and RAM limits. The ledger is
stored in output/, which DPS normally collects after the job finishes.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse

REQUIRED = ("granule", "input_file", "maturity", "date")
OPTIONAL = ("track", "frame", "pass_dir", "matched_stations")
MATURITY = {"provisional", "validated"}
NISAR_FILE = re.compile(r"^NISAR_[LS]2_[A-Z]{2}_GCOV_.*\.(?:h5|hdf5)$", re.I)


def stamp():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def read_manifest(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        missing = set(REQUIRED) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Manifest needs columns {sorted(missing)}")
        result = []
        seen = set()
        for index, row in enumerate(reader, start=2):
            for key in REQUIRED:
                if not (row.get(key) or "").strip():
                    raise ValueError(f"Manifest line {index}: blank {key}")
            maturity = row["maturity"].strip().lower()
            if maturity not in MATURITY:
                raise ValueError(f"Manifest line {index}: only provisional and validated allowed")
            url = row["input_file"].strip()
            name = Path(urlparse(url).path).name
            if not url.startswith("s3://") or not NISAR_FILE.fullmatch(name):
                raise ValueError(f"Manifest line {index}: expected S3 URI to a primary NISAR GCOV .h5: {url}")
            if name.rsplit(".", 1)[0] != row["granule"].strip():
                raise ValueError(f"Manifest line {index}: granule does not match URI basename")
            if row["granule"] in seen:
                raise ValueError(f"Manifest line {index}: duplicate granule {row['granule']}")
            seen.add(row["granule"])
            result.append({**{k: (v or "").strip() for k, v in row.items() if k},
                           "maturity": maturity, "input_file": url})
    if not result:
        raise ValueError("Manifest contains zero granules")
    return result


def dump_atomic(path: Path, rows: list[dict], fieldnames: list[str]):
    temp = path.with_name(path.name + ".tmp")
    with temp.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
        f.flush()
        os.fsync(f.fileno())
    os.replace(temp, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--runner", required=True, type=Path)
    parser.add_argument("--downloader", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("output"))
    parser.add_argument("--base-python", default="python3")
    parser.add_argument("--ease2-res", choices=("50", "100"), default="50")
    parser.add_argument("--backscatter", choices=("sigma0", "gamma0"), default="sigma0")
    parser.add_argument("--mode", choices=("nofilt", "enhanced_lee", "gamma_map", "all"), default="enhanced_lee")
    parser.add_argument("--filter-window", default="7")
    parser.add_argument("--lia", choices=("yes", "no"), default="no")
    parser.add_argument("--frequency", choices=("frequencyA", "frequencyB"), default="frequencyA")
    parser.add_argument("--pol-pair", default="auto")
    parser.add_argument("--workers", default="1")
    parser.add_argument("--tile-workers", default="4")
    parser.add_argument("--max-scenes", type=int, default=0, help="0 = all manifest rows")
    args = parser.parse_args()
    if not args.manifest.is_file():
        parser.error(f"Manifest does not exist: {args.manifest}")
    if not args.runner.is_file() or not args.downloader.is_file():
        parser.error("MAAP runner/downloader is missing")
    # Fail before downloading any data if the manifest isn't complete or well-formed.
    rows = read_manifest(args.manifest)
    if args.max_scenes:
        if args.max_scenes < 0:
            parser.error("max-scenes cannot be negative")
        rows = rows[:args.max_scenes]
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    log_path = out / "single_dps_job_scene_ledger.csv"
    summary_path = out / "single_dps_job_summary.json"
    manifest_copy = out / "input_granule_manifest.csv"
    if args.manifest.resolve() != manifest_copy.resolve():
        manifest_copy.write_bytes(args.manifest.read_bytes())
    fields = ["granule", "maturity", "date", "track", "frame", "pass_dir", "input_file",
              "status", "started_utc", "finished_utc", "error"]
    ledger: list[dict] = []
    settings = {"ease2_res": args.ease2_res, "backscatter": args.backscatter,
                "mode": args.mode, "filter_window": args.filter_window,
                "lia": args.lia, "frequency": args.frequency,
                "clip_to_aoi": "no", "manifest_count": len(rows)}
    print("[SINGLE MAAP JOB] Total scenes:", len(rows), "settings:", settings, flush=True)
    print("[SINGLE MAAP JOB] Full scenes; NO clipping; one GCOV on disk at a time", flush=True)
    with tempfile.TemporaryDirectory(prefix="nisar_gcov_batch_", dir=str(Path.cwd())) as tmp:
        stage = Path(tmp)
        for index, row in enumerate(rows, start=1):
            filename = Path(urlparse(row["input_file"]).path).name
            local = stage / filename
            started = stamp()
            state = "failed"
            detail = ""
            try:
                print(f"[SINGLE MAAP JOB] {index}/{len(rows)}: {filename} {row['maturity']}", flush=True)
                subprocess.run([args.base_python, str(args.downloader), "get",
                                row["input_file"], str(local)], check=True)
                if not local.is_file() or local.stat().st_size == 0:
                    raise RuntimeError("download did not create a nonempty HDF5")
                prior_tiffs = {p.relative_to(out) for p in out.rglob("*.tif")}
                # The existing MAAP runner handles processing, output naming,
                # CRID-specific masking and EASE2 companion creation. It remains
                # foreground and writes every scene to the SAME DPS output folder.
                command = ["bash", str(args.runner), "--input_file", str(local),
                           "--ease2_res", args.ease2_res,
                           "--backscatter", args.backscatter,
                           "--frequency", args.frequency,
                           "--pol_pair", args.pol_pair,
                           "--mode", args.mode,
                           "--filter_window", args.filter_window,
                           "--lia", args.lia,
                           "--workers", args.workers,
                           "--tile_workers", args.tile_workers,
                           "--clip_to_aoi", "no"]
                subprocess.run(command, check=True)
                added_tiffs = {p.relative_to(out) for p in out.rglob("*.tif")} - prior_tiffs
                if not added_tiffs:
                    raise RuntimeError("Processor returned success but created no new GeoTIFF; check logs for skipped input/output-name collision")
                state = "completed"
            except Exception as error:
                detail = f"{type(error).__name__}: {error}"
                print("[SINGLE MAAP JOB] SCENE FAILED:", filename, detail, file=sys.stderr, flush=True)
            finally:
                local.unlink(missing_ok=True)
                ledger.append({**row, "status": state, "started_utc": started,
                               "finished_utc": stamp(), "error": detail})
                dump_atomic(log_path, ledger, fields)
                summary = {"updated_utc": stamp(), "settings": settings,
                           "total": len(rows), "attempted": len(ledger),
                           "completed": sum(x["status"] == "completed" for x in ledger),
                           "failed": sum(x["status"] == "failed" for x in ledger),
                           "remaining": len(rows) - len(ledger),
                           "failures": [x["granule"] for x in ledger if x["status"] == "failed"]}
                summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"[SINGLE MAAP JOB] Finished: {summary['completed']}/{len(rows)} successful; "
          f"{summary['failed']} failures. Check {summary_path}", flush=True)
    # A nonzero exit correctly labels an incomplete MAAP job as failed. MAAP may
    # triage outputs from failed/timed-out jobs, so do not treat it as a guarantee
    # of durable incremental output publication.
    return 0 if summary["failed"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())