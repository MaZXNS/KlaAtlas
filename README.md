# KlaAtlas 数据与分析图

本包使用 2026-11-04 数据库版本 `ROOT-MIN10-CONSERVATIVE-GROUP-20261104`，文献检索截止 2026-10-24，包含 212,444 个位点、568,654 条直接证据、46,280 个蛋白与物种组合及 247 张来源表。SQLite 保存从公开结果表整理的位点、证据、来源、参考序列和绘图输入。

- [数据库与查询说明](database/README.md)
- [分析图](figures/)
- [数据字典](database/DATA_DICTIONARY.tsv)
- [文件校验清单](SHA256SUMS)

## 使用

使用 Python 3.11，绘图字体为 Arial。

```sh
git lfs install
git clone https://github.com/MaZXNS/KlaAtlas.git
cd KlaAtlas
git lfs pull
python verify_package.py
python -m pip install -r requirements.txt
python database/restore_database.py
python verify_package.py --database database/klaatlas.sqlite
python make_figures.py --database database/klaatlas.sqlite --output generated_figures
```

| 分析图 | 内容 |
| --- | --- |
| Figure 1b–f | 图谱覆盖、来源年份、坐标一致性与累积覆盖。 |
| Figure 2 | 位点覆盖、位置分布、聚集、肽段和定位统计。 |
| Figure 3 | 序列特征、对称性与组蛋白分析。 |
| Figure 4 | 亚细胞定位、结构和功能富集。 |

使用 `--figure 1`、`2`、`3`、`4` 选择分析图。`run_metadata.json` 记录输入哈希、运行参数和软件版本。

## 数据

`sites` 与 `evidence` 保存位点及对应来源观察；`proteins`、`sources`、`publications` 和 `sequences` 提供蛋白、来源与序列信息；`figure_data` 保存分析图输入。证据可通过 `source_file_sha256`、`source_table_id`、`source_sheet_or_table` 和 `source_row_number` 追溯。
