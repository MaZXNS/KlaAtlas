#!/usr/bin/env python3
"""Restore the compact KlaAtlas SQLite from its gzip payload."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import tempfile
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gzip', type=Path, default=Path(__file__).resolve().parent / 'klaatlas.sqlite.gz')
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent / 'klaatlas.sqlite')
    parser.add_argument('--sha256', default='ea52a9ebc2ad8e69995f91e353f720e1221562c8c4d001deab911a29db66d171')
    args = parser.parse_args()
    source = args.gzip.resolve()
    target = args.output.resolve()
    if target.exists():
        raise SystemExit('Refusing to overwrite existing database: ' + str(target))
    if not source.is_file():
        raise SystemExit('Missing gzip payload: ' + str(source))
    target.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    handle = tempfile.NamedTemporaryFile(prefix=target.name + '.', suffix='.tmp', dir=target.parent, delete=False)
    temporary = Path(handle.name)
    try:
        with gzip.open(source, 'rb') as reader, handle:
            for chunk in iter(lambda: reader.read(1024 * 1024 * 8), b''):
                digest.update(chunk)
                handle.write(chunk)
        actual = digest.hexdigest()
        if actual != args.sha256:
            raise SystemExit('SHA-256 mismatch: ' + actual)
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()
    print('Restored ' + str(target) + ' (' + str(target.stat().st_size) + ' bytes, sha256 ' + args.sha256 + ')')


if __name__ == '__main__':
    main()
