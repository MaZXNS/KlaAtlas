#!/usr/bin/env python3
"""Generate KlaAtlas analysis figures from the compact database."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import sqlite3
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pymupdf

ROOT = Path(__file__).resolve().parent
CODE = ROOT / 'code'
STATIC = ROOT / 'figures'


def extract_inputs(database: Path, destination: Path, figures: list[str]) -> list[dict[str, str]]:
    """Write the frozen plotting inputs from the database to a temporary tree."""
    inputs = []
    with sqlite3.connect(database.as_uri() + '?mode=ro', uri=True) as connection:
        for figure in figures:
            bucket = f'F{figure}'
            root = destination / bucket
            for name, content in connection.execute(
                    'SELECT name, content FROM figure_data WHERE figure = ? ORDER BY name', (bucket,)):
                target = (root / name).resolve()
                if not target.is_relative_to(root):
                    raise ValueError(f'Unsafe figure input path: {bucket}/{name}')
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(content)
                inputs.append({'input': f'{bucket}/{name}', 'sha256': hashlib.sha256(content).hexdigest()})
    return inputs


def fit_figure1(panel_pdf: Path, output: Path) -> None:
    """Place overview panels on the reference page geometry."""
    with pymupdf.open(panel_pdf) as redrawn, pymupdf.open(STATIC / 'Figure_1.pdf') as reference, pymupdf.open() as combined:
        rect = reference[0].rect
        page = combined.new_page(width=rect.width, height=rect.height)
        page.show_pdf_page(page.rect, redrawn, 0)
        combined.save(output / 'Figure_1.pdf', garbage=4, deflate=True)


def match_figure3_page(figure: Path) -> None:
    """Fit the regenerated vector artwork to the smaller submitted PDF page."""
    with pymupdf.open(STATIC / figure.name) as reference, pymupdf.open(figure) as source, \
            pymupdf.open() as resized:
        rect = reference[0].rect
        page = resized.new_page(width=rect.width, height=rect.height)
        page.show_pdf_page(page.rect, source, 0)
        payload = resized.tobytes(garbage=4, deflate=True)
    figure.write_bytes(payload)


def draw_scientific(data: Path, output: Path, figures: list[str]) -> None:
    sys.path.insert(0, str(CODE))
    import figure1bf
    import figure2
    import figure3
    import figure4

    if '1' in figures:
        figure1bf.render_panels_bf(data / 'F1', data)
        fit_figure1(data / 'Figure_1b_f.pdf', output)
    if '2' in figures:
        figure2.draw(data / 'F2', output)
    if '3' in figures:
        figure3.draw(data / 'F3', output)
        match_figure3_page(output / 'Figure_3.pdf')
    if '4' in figures:
        figure4.draw(data / 'F4', output)


def main() -> None:
    """Draw selected figures without replacing the frozen source assets."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--database', type=Path, default=ROOT / 'database' / 'klaatlas.sqlite')
    parser.add_argument('--output', type=Path, default=ROOT / 'generated_figures')
    parser.add_argument('--figure', choices=['all', '1', '2', '3', '4'], default='all')
    args = parser.parse_args()

    database = args.database.resolve()
    output = args.output.resolve()
    if output == STATIC or output.is_relative_to(STATIC):
        parser.error('Use an output directory outside the frozen figures directory')
    output.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR', str(output / '.mplconfig'))
    os.environ.setdefault('MPLBACKEND', 'Agg')
    figures = ['1', '2', '3', '4'] if args.figure == 'all' else [args.figure]

    inputs = []
    if figures:
        with tempfile.TemporaryDirectory(prefix='klaatlas_plot_', dir=output) as directory:
            data = Path(directory)
            inputs = extract_inputs(database, data, figures)
            draw_scientific(data, output, figures)

    receipt = {
        'run_utc': datetime.now(timezone.utc).isoformat(),
        'database': str(database) if figures else None,
        'output': str(output),
        'figure': args.figure,
        'redrawn': figures,
        'inputs': inputs,
        'python': sys.version.split()[0],
        'packages': {name: importlib.metadata.version(name) for name in
                     ['numpy', 'pandas', 'matplotlib', 'logomaker', 'PyMuPDF', 'Pillow']},
    }
    (output / 'run_metadata.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print('Wrote figures to ' + str(output))


if __name__ == '__main__':
    main()
