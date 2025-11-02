"""
Example script to run triaxial stress-strain analysis

This script demonstrates how to use the TriaxialAnalyzer class to:
1. Load and process triaxial compression test data
2. Generate summary statistics
3. Create visualizations
4. Save results
"""

from triaxial_analyzer import TriaxialAnalyzer
import sys
from pathlib import Path

def main():
    """Main analysis workflow"""
    
    # Configuration
    DATA_FILE = '../data/your_data.csv'  # Update this path
    OUTPUT_DIR = Path('../output')
    PLOTS_DIR = OUTPUT_DIR / 'plots'
    STATS_DIR = OUTPUT_DIR / 'statistics'
    
    # Analysis parameters
    MAX_STRESS_FRACTION = 0.70  # Analyze up to 70% of max stress
    SMOOTH_WINDOW = 11          # Smoothing window (must be odd)
    SMOOTH_ORDER = 2            # Smoothing polynomial order
    
    print("="*80)
    print("TRIAXIAL STRESS-STRAIN ANALYSIS")
    print("="*80)
    print(f"Data file: {DATA_FILE}")
    print(f"Analysis range: 0-{MAX_STRESS_FRACTION*100:.0f}% of max stress")
    print(f"Smoothing: window={SMOOTH_WINDOW}, order={SMOOTH_ORDER}")
    print("="*80)
    
    # Initialize analyzer
    try:
        analyzer = TriaxialAnalyzer(
            filepath=DATA_FILE,
            max_stress_fraction=MAX_STRESS_FRACTION,
            smooth_window=SMOOTH_WINDOW,
            smooth_order=SMOOTH_ORDER
        )
    except Exception as e:
        print(f"\nERROR: Could not load data file: {e}")
        print(f"Please check that '{DATA_FILE}' exists and is formatted correctly.")
        sys.exit(1)
    
    # Process all samples
    print("\n" + "="*80)
    print("PROCESSING SAMPLES")
    print("="*80)
    analyzer.process_all_samples()
    
    # Generate and save summary statistics
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    summary = analyzer.get_summary_statistics()
    print(summary.to_string(index=False))
    
    # Save summary to CSV
    summary_file = STATS_DIR / 'summary_statistics.csv'
    summary.to_csv(summary_file, index=False)
    print(f"\n✓ Saved summary statistics to: {summary_file}")
    
    # Generate plots for all samples
    print("\n" + "="*80)
    print("GENERATING PLOTS")
    print("="*80)
    
    sample_ids = list(analyzer.results.keys())
    
    # Individual sample plots
    print(f"\nCreating individual sample plots ({len(sample_ids)} samples)...")
    for i, sample_id in enumerate(sample_ids, 1):
        # Create safe filename (replace special characters)
        safe_filename = sample_id.replace('/', '-').replace('\\', '-')
        plot_file = PLOTS_DIR / f'sample_{safe_filename}.png'
        
        print(f"  [{i}/{len(sample_ids)}] {sample_id}")
        analyzer.plot_sample(sample_id, save_path=str(plot_file))
    
    # Multi-sample comparison plot
    print("\nCreating multi-sample comparison plot...")
    comparison_file = PLOTS_DIR / 'comparison_all_samples.png'
    analyzer.plot_comparison(save_path=str(comparison_file))
    print(f"✓ Saved to: {comparison_file}")
    
    # Statistical summary plot
    print("\nCreating statistical summary plot...")
    stats_plot_file = PLOTS_DIR / 'statistics_summary.png'
    analyzer.plot_statistics_summary(save_path=str(stats_plot_file))
    print(f"✓ Saved to: {stats_plot_file}")
    
    # Summary
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"Processed: {len(sample_ids)} samples")
    print(f"Plots saved to: {PLOTS_DIR}")
    print(f"Statistics saved to: {STATS_DIR}")
    print("="*80)

if __name__ == "__main__":
    main()