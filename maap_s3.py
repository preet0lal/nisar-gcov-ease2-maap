#!/usr/bin/env python3
"""Temporary NISAR (ASF) Earthdata S3 access for MAAP DPS jobs.

Run with the MAAP base-image Python (it has maap-py and s3fs), not the
nisar_gcov_ease2 conda environment.

  maap_s3.py env            -> prints shell `export` lines for AWS_* variables
                               (eval them; never echo them to logs)
  maap_s3.py get URL DEST   -> downloads one s3:// object to DEST
"""
from __future__ import annotations

import shlex
import sys
from pathlib import Path

S3_CREDENTIAL_URL = "https://nisar.asf.earthdatacloud.nasa.gov/s3credentials"


def credentials() -> dict[str, str]:
    from maap.maap import MAAP

    creds = MAAP().aws.earthdata_s3_credentials(S3_CREDENTIAL_URL)
    required = ("accessKeyId", "secretAccessKey", "sessionToken")
    if not isinstance(creds, dict) or any(not creds.get(k) for k in required):
        raise SystemExit("ERROR: MAAP did not return usable Earthdata S3 credentials")
    return {
        "AWS_ACCESS_KEY_ID": creds["accessKeyId"],
        "AWS_SECRET_ACCESS_KEY": creds["secretAccessKey"],
        "AWS_SESSION_TOKEN": creds["sessionToken"],
        "AWS_REGION": "us-west-2",
        "AWS_DEFAULT_REGION": "us-west-2",
        "AWS_REQUEST_PAYER": "requester",
    }


def download(url: str, destination: str) -> None:
    import s3fs

    env = credentials()
    fs = s3fs.S3FileSystem(
        key=env["AWS_ACCESS_KEY_ID"], secret=env["AWS_SECRET_ACCESS_KEY"],
        token=env["AWS_SESSION_TOKEN"], requester_pays=True,
        client_kwargs={"region_name": "us-west-2"},
    )
    target = Path(destination)
    target.parent.mkdir(parents=True, exist_ok=True)
    print(f"[MAAP] Downloading {url} -> {target}", file=sys.stderr, flush=True)
    fs.get(url, str(target))
    if not target.is_file() or target.stat().st_size == 0:
        raise SystemExit(f"ERROR: download empty or missing: {target}")
    print(f"[MAAP] Downloaded {target.stat().st_size / 1024**3:.2f} GiB",
          file=sys.stderr, flush=True)


def main(argv: list[str]) -> None:
    if argv[:1] == ["env"] and len(argv) == 1:
        for key, value in credentials().items():
            print(f"export {key}={shlex.quote(value)}")
    elif argv[:1] == ["get"] and len(argv) == 3:
        download(argv[1], argv[2])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
