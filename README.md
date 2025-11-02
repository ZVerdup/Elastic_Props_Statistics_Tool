# Triaxial Stress-Strain Analysis Tool

A Python-based tool for analyzing stress-strain data from rock triaxial compression tests. Computes continuous profiles of Young's Modulus and Poisson's Ratio and enables statistical comparisons between samples.

## Features

- **Multi-sample processing**: Analyze multiple samples from a single data file
- **Flexible analysis range**: Configure the stress threshold for analysis (default: 70% of max deviatoric stress)
- **Orientation-aware**: Handles both Horizontal and Vertical sample orientations
- **Comprehensive metrics**:
  - Young's Modulus (YM) profiles and statistics
  - Poisson's Ratio A (PRA) from radial strain A
  - Poisson's Ratio B (PRB) from radial strain B
  - Average Poisson's Ratio (PRAvg) for Vertical samples
- **Data smoothing**: Savitzky-Golay filtering to handle noisy derivatives
- **Rich visualizations**:
  - Individual sample plots (stress-strain curves + moduli evolution)
  - Multi-sample comparison plots
  - Statistical summary box plots

## Project Structure

```
triaxial-stress-analysis/
├── src/
│   ├── triaxial_analyzer.py    # Main analysis class
│   └── run_analysis.py          # Example script to run analysis
├── data/
│   └── (place your CSV files here)
├── output/
│   ├── plots/                   # Generated plots saved here
│   └── statistics/              # Summary statistics CSV files
├── tests/
│   └── (unit tests)
├── docs/
│   └── (documentation)
├── requirements.txt
├── README.md
└── .gitignore
```

## Installation

1. **Clone the repository** (or create it locally):
```bash
git clone <your-repo-url>
cd triaxial-stress-analysis
```

2. **Install dependencies**:
```bash
pip install -r requirements.txt
```

## Usage

### Basic Analysis

```python
from src.triaxial_analyzer import TriaxialAnalyzer

# Initialize analyzer
analyzer = TriaxialAnalyzer(
    filepath='data/your_data.csv',
    max_stress_fraction=0.70,  # Analyze up to 70% of max stress
    smooth_window=11,          # Smoothing window size
    smooth_order=2             # Polynomial order for smoothing
)

# Process all samples
analyzer.process_all_samples()

# Get summary statistics
summary = analyzer.get_summary_statistics()
print(summary)

# Save summary
summary.to_csv('output/statistics/summary.csv', index=False)
```

### Plotting Individual Samples

```python
# Plot a specific sample
analyzer.plot_sample('1-107-1H', save_path='output/plots/sample_1-107-1H.png')
```

### Multi-Sample Comparison

```python
# Compare all samples
analyzer.plot_comparison(save_path='output/plots/comparison_all.png')

# Compare specific samples
analyzer.plot_comparison(
    sample_ids=['1-107-1H', '2-108-1H', '3-109-1H'],
    save_path='output/plots/comparison_selected.png'
)
```

### Statistical Summary

```python
# Create box plots of statistics
analyzer.plot_statistics_summary(save_path='output/plots/statistics_summary.png')
```

## Data Format

The tool expects CSV files with the following columns:

- `sample_id`: Unique identifier for each sample
- `orientation`: "Horizontal" or "Vertical"
- `sigma_d_psi`: Deviatoric stress (psi)
- `eps_ax`: Axial strain
- `eps_rad_a`: Radial strain A
- `eps_rad_b`: Radial strain B

Additional columns are preserved but not required for analysis.

## Configuration

### Analysis Parameters

- **max_stress_fraction**: Fraction of maximum deviatoric stress to include in analysis (default: 0.70)
- **smooth_window**: Window size for Savitzky-Golay smoothing filter (default: 11, must be odd)
- **smooth_order**: Polynomial order for smoothing (default: 2)

### Adjusting Parameters

```python
# For noisier data, increase smoothing
analyzer = TriaxialAnalyzer(
    filepath='data/noisy_data.csv',
    smooth_window=15,
    smooth_order=3
)

# To analyze a different stress range
analyzer = TriaxialAnalyzer(
    filepath='data/your_data.csv',
    max_stress_fraction=0.50  # Analyze up to 50% of max stress
)
```

## Output

### Summary Statistics CSV

Includes for each sample:
- Sample ID and orientation
- Maximum stress and number of analysis points
- Young's Modulus: mean, std, median, min, max
- Poisson's Ratios (A, B, and Avg for Vertical): mean, std, median

### Plots

1. **Individual sample plots**: 2x2 grid showing:
   - Stress-strain curve with analysis range
   - Young's Modulus evolution
   - Poisson's Ratio A evolution
   - Poisson's Ratio B (or Avg for Vertical) evolution

2. **Comparison plots**: Multi-sample overlay of all profiles

3. **Statistics plots**: Box plots summarizing distributions across samples

## Contributing

Contributions welcome! Please feel free to submit pull requests or open issues.

## License

[Specify your license here]

## Contact

[Your contact information]

## Acknowledgments

Developed for analysis of rock mechanics triaxial compression test data.