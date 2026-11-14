# 修回稿便携数据库与图件复现

2026 年 11 月 5 日修回稿使用 212,444 个位点、568,654 条直接证据、46,280 个蛋白与物种组合、144 份来源报告和 247 张来源表。该包的 SQLite 提供位点、证据、来源和冻结绘图输入。本包固定使用 2026-11-04 数据库版本 `ROOT-MIN10-CONSERVATIVE-GROUP-20261104`，文献检索截止 2026-10-24。所含数据是从作者公开结果表整理、标准化后的客观证据提取，不是原始质谱数据。

- [数据库、常用查询与来源](database/README.md)
- [定稿 PDF 副本](figures/)
- [本目录文件校验](SHA256SUMS)

## 恢复与重绘

使用 Python 3.11 和 `requirements.txt`，系统应能加载 Arial 字体。解压得到 1,910,112,256 字节的 SQLite。可用 `--output` 将它恢复到临时位置，并把同一路径传给 `make_figures.py` 的 `--database`。

```sh
git lfs install
git clone https://github.com/MaZXNS/KlaAtlas.git
cd KlaAtlas
git lfs pull
python verify_package.py
python -m pip install -r requirements.txt
python database/restore_database.py
python verify_package.py --database database/klaatlas.sqlite
python make_figures.py --database database/klaatlas.sqlite --output /path/to/output
```

| 图件 | 复现方式 |
| --- | --- |
| Figure 1 | 原样嵌入作者定稿面板 a。从冻结统计输入重绘 b–f，保留定稿的两页尺寸。 |
| Figure 2 | 从冻结的覆盖、位置、聚集、肽段和定位统计表重绘。 |
| Figure 3 | 从冻结的序列、对称性和组蛋白结果重绘，PDF 匹配投稿页面尺寸。 |
| Figure 4 | 从冻结的定位、结构和富集结果重绘。 |
| Figure 5、Supporting Information、TOC | 复制最终 PDF，Figure 5 使用投稿中的网站截图版。 |

Figure 1a 定稿显示 7,531 条去重记录、6,068 条初筛排除、1,463 条详审、1,319 条详审排除、144 份来源报告和 247 张整合来源表。面板 a 的来源和哈希见 `metadata.json`。`--figure` 支持 `1` 至 `5` 和 `all`，输出的 `run_metadata.json` 记录输入、运行参数及软件版本。

## 数据库范围

`figure_data` 保存 27 项冻结绘图输入。便携库保留 7,139 条文献筛选记录。`final_publication_class='included_or_supplementary_tier'` 返回纳入或补充层记录，144 是修回稿的来源报告口径。证据表保留 `source_file_sha256`、来源文件、来源表和原始行号。模型拟合、结构、结构域与功能注释使用科研母库及长期归档中的完整统计包。

主图 1a 保留定稿流程图；数据驱动的主图 1b–f、2–4 从冻结输入重绘。Figure 5、含补充图的完整 Supporting Information 和 TOC 提供定稿 PDF。此最小包不重新执行统计模型拟合。
