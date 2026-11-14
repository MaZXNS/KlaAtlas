# Database

```text
klaatlas.sqlite.gz   Gzip-compressed SQLite database (about 342 MiB; 1.91 GB unpacked)
restore_database.py  Restore the database and verify its SHA-256 checksum
```

From the package root:

```sh
python database/restore_database.py
```

The default output is `database/klaatlas.sqlite`, with SHA-256:

```text
ea52a9ebc2ad8e69995f91e353f720e1221562c8c4d001deab911a29db66d171
```

## Tables

| Table | Contents |
| --- | --- |
| `sites` | 212,444 rows, one per Kla site. |
| `evidence` | 568,654 rows, one per source observation. |
| `proteins` | 46,280 protein–taxon entries. |
| `publications` | 7,139 screening records. Use `final_publication_class = 'included_or_supplementary_tier'` to select included or supplementary records. |
| `sources` | 247 active source tables. |
| `sequences` | 34,485 reference sequences. |
| `figure_data` | 26 frozen analysis figure inputs stored as BLOBs. |

The [data dictionary](DATA_DICTIONARY.tsv) lists the database schema. `evidence.source_file_sha256` identifies source file content; `source_table_id`, `source_sheet_or_table` and `source_row_number` locate the corresponding source records.

## Example queries

From the package root:

```sh
sqlite3 -readonly database/klaatlas.sqlite 'SELECT COUNT(*) FROM sites; SELECT COUNT(*) FROM evidence;'
sqlite3 -readonly database/klaatlas.sqlite 'PRAGMA table_info(evidence);'
```

For interactive exploration, visit [KlaAtlas](https://www.kla-atlas.com/).
