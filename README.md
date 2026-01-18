# Enzyme Correlator

[![CI](https://github.com/JPBureik/enzyme_correlator/actions/workflows/ci.yml/badge.svg)](https://github.com/JPBureik/enzyme_correlator/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JPBureik/enzyme_correlator/branch/master/graph/badge.svg)](https://codecov.io/gh/JPBureik/enzyme_correlator)
[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Checked with mypy](https://www.mypy-lang.org/static/mypy_badge.svg)](https://mypy-lang.org/)

Correlation matrix calculator for enzyme activity analysis. This tool calculates the correlation matrix for enzyme activity with respect to substrates and plots the result for graphical quantitative analysis.

This software was developed as supplementary material for [Sharma et al. (2021)](https://pubmed.ncbi.nlm.nih.gov/33909340/).

## Installation

```bash
git clone https://github.com/JPBureik/enzyme_correlator.git
cd enzyme_correlator
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install -e .
```

For development:

```bash
pip install -e ".[dev]"
pre-commit install
```

## Usage

### Running the Application

After installation:

```bash
enzyme-correlator
```

Or using Python directly:

```bash
python enzyme_correlations.py
```

### Data Format

The input data must be a CSV file with the following format:

- Semicolon-delimited (`;`)
- First row: substrate/condition names (header)
- First column: enzyme names
- Data cells: decimal numbers using comma as decimal separator (e.g., `0,85`)

Example structure:

```csv
;Substrate1;Substrate2;Substrate3
Enzyme1;0,90;0,85;0,80
Enzyme2;0,88;0,92;0,78
Enzyme3;0,20;0,25;0,30
```

To create a correctly formatted CSV file from an Excel file:

1. Delete all superfluous rows and columns (including e.g., substrate pictograms)
2. Save as CSV with enzymes as rows and substrates as columns
3. Each enzyme and substrate should have one header cell containing its name

### Application Features

- **Load Data**: Import CSV files with enzyme activity data
- **Show Enzyme Grouping**: Display enzymes grouped by correlation
- **Plot Correlation Matrix**: Visualize correlations as a heatmap
- **Plot Histogram**: Show distribution of correlation values
- **Grouping Cutoff Slider**: Adjust the correlation threshold for grouping (default: 0.85)
- **Save Figure**: Export visualizations to image files

## Development

### Running Tests

```bash
pytest tests/ -v
```

### Running Linters

```bash
ruff check src/
mypy src/
```

### Pre-commit Hooks

Pre-commit hooks are configured for code quality:

```bash
pre-commit run --all-files
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Citation

If you use this software in your research, please cite:

> Sharma et al. (2021). Conversion of five proluciferin esters by human cytochrome P450 enzymes. *Biotechnology Journal*. https://pubmed.ncbi.nlm.nih.gov/33909340/
