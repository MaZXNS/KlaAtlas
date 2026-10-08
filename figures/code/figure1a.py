from __future__ import annotations



import json


import textwrap

from pathlib import Path


import matplotlib.pyplot as plt

from matplotlib.patches import Rectangle

from plotting_common import export


LINE = "#7fa4b6"

FILL = "#e7f1f6"

TEXT = "#1a1a1a"

MUTED = "#617783"

TEAL = "#258975"

PURPLE = "#7863a6"

SIDE_LABEL_SIZE_PT = 8.7

SIDE_LABEL_WEIGHT = "normal"

def build(args: object) -> None:
    """Render the isolated panel and write numerical and geometry QA records."""
    args.output.mkdir(parents=True, exist_ok=True)

    raw = json.loads(args.source_text.read_text())
    assert len(raw) == 12
    raw[2] = "Records after deduplication and identity reconciliation\n(n = 7,139)"
    raw[11] = "INCLUSION"

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 10,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "axes.linewidth": 0.6,
        }
    )
    fig = plt.figure(figsize=(178 / 25.4, 184 / 25.4), facecolor="white")
    # Preserve the stored physical size across Matplotlib figure-manager rounding.
    fig.set_size_inches(178 / 25.4, 184 / 25.4, forward=False)
    axis = fig.add_axes([0, 0, 1, 1])
    axis.set_xlim(0, 178)
    axis.set_ylim(184, 0)
    axis.axis("off")
    containment_items: list[tuple[Rectangle, list[object]]] = []

    def text(
        x: float,
        y: float,
        value: str,
        size: float = 10,
        horizontal: str = "center",
        weight: str = "normal",
        color: str = TEXT,
        rotation: int = 0,
    ) -> object:
        return axis.text(
            x,
            y,
            value,
            fontsize=size,
            ha=horizontal,
            va="center",
            fontweight=weight,
            color=color,
            linespacing=1.2,
            rotation=rotation,
            rotation_mode="anchor",
        )

    def box(
        x: float,
        y: float,
        width: float,
        height: float,
        body: str,
        fill: str = FILL,
        size: float = 10,
        bold: bool = False,
    ) -> None:
        patch = Rectangle(
            (x, y),
            width,
            height,
            facecolor=fill,
            edgecolor=TEAL if bold else LINE,
            lw=0.8,
        )
        axis.add_patch(patch)
        item = text(
            x + width / 2,
            y + height / 2,
            body,
            size,
            weight="bold" if bold else "normal",
            color=TEAL if bold else TEXT,
        )
        containment_items.append((patch, [item]))

    def reasons(x: float, y: float, width: float, height: float, body: str) -> None:
        parts = body.split("\n")
        header = parts[0].replace(" (n =", "\n(n =")
        patch = Rectangle(
            (x, y), width, height, facecolor="#f6f3fa", edgecolor="#b1a5c5", lw=0.7
        )
        axis.add_patch(patch)
        lines: list[str] = []
        for line in parts[1:]:
            if line.strip().startswith("Processed site-level table unavailable:"):
                lines += ["Processed site-level table", "unavailable: 11"]
            else:
                lines += textwrap.wrap(
                    line.strip(), width=40, break_long_words=False, break_on_hyphens=False
                )
        line_mm = 9.4 * 1.2 * 25.4 / 72
        total_height = 2 * line_mm + 2.3 + len(lines) * line_mm
        top = y + (height - total_height) / 2
        items = [text(x + 3, top + line_mm, header, 9.4, "left", "bold", PURPLE)]
        items.append(
            text(
                x + 3,
                top + 2 * line_mm + 2.3 + len(lines) * line_mm / 2,
                "\n".join(lines),
                9.4,
                "left",
            )
        )
        containment_items.append((patch, items))

    def arrow(x: float, y: float, xx: float, yy: float) -> None:
        axis.annotate(
            "",
            xy=(xx, yy),
            xytext=(x, y),
            arrowprops={
                "arrowstyle": "-|>",
                "color": MUTED,
                "lw": 0.65,
                "mutation_scale": 6,
                "shrinkA": 0,
                "shrinkB": 0,
            },
        )

    text(3, 5, raw[0], 12, "left", "bold", "#000000")
    box(19, 8, 155, 26, raw[1], size=9.6)
    box(
        19,
        43,
        73,
        21,
        "Records after deduplication\nand identity reconciliation\n(n = 7,139)",
        size=10.5,
    )
    box(19, 77, 73, 20, raw[3].replace("and metadata", "and\nmetadata"), size=11)
    reasons(101, 67.5, 73, 39, raw[4].replace("1,846No", "1,846\nNo"))
    box(
        19,
        122,
        73,
        20,
        raw[5].replace("evidence assessment", "evidence\nassessment"),
        size=11,
    )
    reasons(101, 110.5, 73, 43, raw[6])
    box(
        19,
        162,
        155,
        18,
        raw[7].replace("，", ","),
        fill="#e7f3ef",
        size=10.7,
        bold=True,
    )
    arrow(55.5, 34, 55.5, 43)
    arrow(55.5, 64, 55.5, 77)
    arrow(92, 87, 101, 87)
    arrow(55.5, 97, 55.5, 122)
    arrow(92, 132, 101, 132)
    arrow(55.5, 142, 55.5, 162)

    # Every side label uses the same typography. Each guide-line midpoint equals
    # the center of the corresponding module group. The inclusion guide is
    # lengthened symmetrically around the retained-module center (171 mm).
    side_labels = [
        {"text": "IDENTIFICATION", "line_start": 8.0, "line_end": 64.0, "module_center": 36.0},
        {"text": "SCREENING", "line_start": 67.5, "line_end": 106.5, "module_center": 87.0},
        {"text": "EVIDENCE ASSESSMENT", "line_start": 110.5, "line_end": 153.5, "module_center": 132.0},
        {"text": "INCLUSION", "line_start": 158.5, "line_end": 183.5, "module_center": 171.0},
    ]
    label_artists = []
    for entry in side_labels:
        midpoint = (entry["line_start"] + entry["line_end"]) / 2
        assert abs(midpoint - entry["module_center"]) < 1e-12
        axis.plot([14, 14], [entry["line_start"], entry["line_end"]], color=LINE, lw=0.8)
        artist = text(
            9,
            midpoint,
            entry["text"],
            SIDE_LABEL_SIZE_PT,
            weight=SIDE_LABEL_WEIGHT,
            color=MUTED,
            rotation=90,
        )
        label_artists.append(artist)

    export(fig, args.output / "Figure_1a")
    plt.close(fig)

def draw(data: Path, output: Path) -> None:
    """Draw the fixed selection workflow."""
    from types import SimpleNamespace
    build(SimpleNamespace(source_text=data/"selection.json",output=output))
