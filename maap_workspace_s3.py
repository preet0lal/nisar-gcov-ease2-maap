#!/usr/bin/env python3
"""Fetch a user-owned MAAP workspace S3 object with temporary workspace credentials.

Used only when DPS did not stage a File Input manifest automatically.
Never print credentials to stdout or logs.
"""
from __future__ import annotations

import sys
from pathlib import Path
from urllib.parse import urlparse


def main(argv):
    if len(argv) != 3 or argv[0] != 'get':
        raise SystemExit('Usage: maap_workspace_s3.py get s3://bucket/path/file.csv DEST')
    uri, dest = argv[1:]
    u = urlparse(uri)
    if u.scheme != 's3' or not u.netloc or not u.path.lstrip('/'):
        raise SystemExit('Expected valid s3:// workspace URI')
    from maap.maap import MAAP
    import s3fs
    reply = MAAP().aws.workspace_bucket_credentials()
    credentials = reply['credentials']
    paths = reply.get('authorized_s3_paths', [])
    key = u.path.lstrip('/')
    accessible = any(
        p.get('bucket') == u.netloc and
        (key == p.get('prefix','').strip('/') or
         key.startswith(p.get('prefix','').strip('/') + '/') or
         not p.get('prefix','').strip('/'))
        for p in paths
    )
    if not accessible:
        raise SystemExit(f'Workspace credentials do not authorize {uri}')
    fs = s3fs.S3FileSystem(
        key=credentials['aws_access_key_id'],
        secret=credentials['aws_secret_access_key'],
        token=credentials['aws_session_token'],
        client_kwargs={'region_name':'us-west-2'},
    )
    target = Path(dest)
    target.parent.mkdir(parents=True, exist_ok=True)
    fs.get(f'{u.netloc}/{key}', str(target))
    if not target.is_file() or target.stat().st_size == 0:
        raise SystemExit(f'Staged workspace input missing or empty: {target}')
    print(f'[MAAP] Staged workspace manifest {target.name}', file=sys.stderr)


if __name__ == '__main__':
    main(sys.argv[1:])