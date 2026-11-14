# database

~~~text
klaatlas.sqlite.gz   1.91 GB SQLite 的 gzip（345 MB）
restore_database.py  解压并校验 SHA-256
~~~

恢复：

~~~sh
python restore_database.py
~~~

默认输出 database/klaatlas.sqlite，SHA-256：
06f229562db8b9da9f15d1a9281d8227b18de99442220d4987a186999c3b1b51

恢复后可用 sqlite3 查询。常用表：

- sites：212,444 行，每个 Kla 位点一行。
- evidence：568,654 行，每条来源观察一行。
- proteins：46,280 行，蛋白 × 物种汇总。
- publications：7,139 行筛选记录；纳入记录过滤 final_publication_class = 'included_or_supplementary_tier'。
- sources：247 行当前活跃来源表。
- sequences：34,485 条参考序列。
- figure_data：Figure 1-4 的冻结绘图输入（BLOB）。

[数据字典](DATA_DICTIONARY.tsv)列出实际表结构。`evidence.source_file_sha256` 连到项目来源对象清单，`source_table_id`、`source_sheet_or_table`、`source_row_number` 可用于定位对应原始记录。模型、结构、结构域或功能注释数据请从科研母库或长期归档获取。

~~~sh
sqlite3 -readonly klaatlas.sqlite 'SELECT COUNT(*) FROM sites; SELECT COUNT(*) FROM evidence;'
sqlite3 -readonly klaatlas.sqlite 'PRAGMA table_info(evidence);'
~~~
