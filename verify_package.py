#!/usr/bin/env python3
"""Verify package checksums and optionally the restored scientific extract."""
from __future__ import annotations
import argparse
import hashlib
import json
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    """Hash a file in bounded memory."""
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def verify(root: Path, database: Path | None = None) -> dict:
    """Check every listed payload; fail on corrupt files or wrong baseline."""
    root = root.resolve()
    metadata = json.loads((root / 'metadata.json').read_text())
    checked = []
    for line in (root / 'SHA256SUMS').read_text().splitlines():
        expected, name = line.split('  ', 1)
        path = (root / name).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise ValueError(f'Missing or unsafe file: {name}')
        if sha256(path) != expected:
            raise ValueError(f'Checksum mismatch: {name}')
        checked.append(name)
    result = {'status': 'PASS', 'checked_files': len(checked),
              'source_database': metadata['source_database']}
    if database is not None:
        database = database.resolve()
        if sha256(database) != metadata['minimal_database']['sqlite_sha256']:
            raise ValueError('Restored database does not match the packaged extract')
        with sqlite3.connect(database.as_uri() + '?mode=ro', uri=True) as con:
            if con.execute('PRAGMA quick_check').fetchone()[0] != 'ok':
                raise ValueError('SQLite quick_check failed')
            counts = {}
            for table, spec in metadata['minimal_database']['tables'].items():
                if not table.replace('_', '').isalnum():
                    raise ValueError('Unsafe table name')
                counts[table] = con.execute(f'SELECT COUNT(*) FROM "{table}"').fetchone()[0]
                if counts[table] != spec['rows']:
                    raise ValueError(f'Row count mismatch: {table}')
            result['database_counts'] = counts
    return result


def main() -> None:
    """Validate a package without writing or changing scientific data."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--database', type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.root, args.database), indent=2))


if __name__ == '__main__':
    main()
