"""Load the selected analysis inputs from SQLite."""
from __future__ import annotations
import csv,gzip,hashlib,sqlite3,zlib
from pathlib import Path


def restore(database: Path, destination: Path, paths: list[str]|None=None) -> list[dict]:
    """Read the registered inputs and core-derived records."""
    manifest=[]
    with sqlite3.connect(f'file:{database.resolve()}?mode=ro',uri=True) as connection:
        for path,digest,codec,payload in connection.execute('SELECT f.original_path,f.sha256,c.storage_codec,c.payload FROM analysis_files f JOIN frozen_content c ON c.content_sha256=f.sha256 ORDER BY f.original_path'):
            if paths is not None and path not in paths:continue
            relative=Path(path)
            if relative.is_absolute() or '..' in relative.parts:raise ValueError('Invalid input registry path')
            content=zlib.decompress(payload) if codec=='zlib' else payload
            if hashlib.sha256(content).hexdigest()!=digest:raise ValueError(f'Input content changed: {path}')
            target=destination/relative;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(content)
            manifest.append({'input':path,'source_sha256':digest,'representation':'frozen_exact_content'})
        for path,query,columns,row_count,digest in connection.execute('SELECT input_path,sql_query,columns_json,row_count,source_content_sha256 FROM analysis_exports'):
            if paths is not None and path not in paths:continue
            target=destination/path;target.parent.mkdir(parents=True,exist_ok=True);cursor=connection.execute(query)
            with gzip.open(target,'wt',newline='') as handle:
                writer=csv.writer(handle,delimiter='\t',lineterminator='\n');writer.writerow([x[0] for x in cursor.description]);writer.writerows(cursor)
            manifest.append({'input':path,'source_sha256':digest,'representation':'generated_from_normalized_core','rows':row_count})
    if paths is not None:
        missing=set(paths)-{row['input'] for row in manifest}
        if missing:raise FileNotFoundError('Unregistered scientific inputs: '+', '.join(sorted(missing)))
    return manifest
