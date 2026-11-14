# 修回稿绘图模块

从包根目录运行 `python make_figures.py --output /path/to/output`。`figure1bf.py` 与 `coverage_panel.py` 来自修回稿实际绘图的冻结代码快照，重绘 Figure 1 面板 b–f。覆盖模块沿用获核验的坐标和数值，按最终投稿的 Figure 1c 使用 log10 来源报告轴、0.5 起点及 1、10、100 刻度。`make_figures.py` 拼接作者定稿面板 a，调用 `figure2.py` 至 `figure4.py` 重绘其余数据图，并复制最终 Figure 5、SI 和 TOC。`plotting_common.py` 与 `figure_alignment.py` 校验导出几何，`functional_panels.py` 与 `heatmap_style.py` 绘制 Figure 4 面板。
