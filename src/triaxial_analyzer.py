"""
Triaxial Compression Test Analysis Tool
Analyzes stress-strain data from rock triaxial compression tests
Computes Young's Modulus and Poisson's Ratio profiles
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import savgol_filter
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

class TriaxialAnalyzer:
    """
    Analyzes triaxial compression test data for multiple samples
    """
    
    def __init__(self, filepath, max_stress_fraction=0.70, smooth_window=11, smooth_order=2):
        """
        Initialize analyzer
        
        Parameters:
        -----------
        filepath : str
            Path to CSV file with triaxial data
        max_stress_fraction : float
            Fraction of max deviatoric stress to analyze (default 0.70 = 70%)
        smooth_window : int
            Window size for Savitzky-Golay smoothing (must be odd)
        smooth_order : int
            Polynomial order for smoothing
        """
        self.filepath = filepath
        self.max_stress_fraction = max_stress_fraction
        self.smooth_window = smooth_window
        self.smooth_order = smooth_order
        
        # Load and parse data
        self.df = self._load_data()
        self.samples = {}
        self.results = {}
        
        print(f"Loaded {len(self.df)} rows")
        print(f"Found {self.df['sample_id'].nunique()} unique samples")
        
    def _load_data(self):
        """Load CSV data, trying different delimiters"""
        try:
            # Try comma first
            df = pd.read_csv(self.filepath)
            if len(df.columns) == 1:
                # Try tab
                df = pd.read_csv(self.filepath, sep='\t')
            if len(df.columns) == 1:
                # Try space
                df = pd.read_csv(self.filepath, sep=r'\s+')
            return df
        except Exception as e:
            raise ValueError(f"Could not load file: {e}")
    
    def _smooth_data(self, data):
        """Apply Savitzky-Golay filter to smooth noisy data"""
        if len(data) < self.smooth_window:
            return data
        try:
            return savgol_filter(data, self.smooth_window, self.smooth_order)
        except:
            return data
    
    def _compute_derivative(self, y, x):
        """Compute derivative dy/dx using central differences"""
        # Remove any NaN values
        mask = ~(np.isnan(x) | np.isnan(y))
        x_clean = x[mask]
        y_clean = y[mask]
        
        if len(x_clean) < 2:
            return np.full_like(x, np.nan)
        
        # Compute derivative using numpy gradient
        dydx = np.gradient(y_clean, x_clean)
        
        # Map back to original indices
        result = np.full_like(x, np.nan, dtype=float)
        result[mask] = dydx
        
        return result
    
    def process_sample(self, sample_id):
        """
        Process a single sample: extract data, compute moduli, analyze
        
        Parameters:
        -----------
        sample_id : str
            Sample identifier
            
        Returns:
        --------
        dict : Analysis results for this sample
        """
        # Extract sample data
        sample_data = self.df[self.df['sample_id'] == sample_id].copy()
        sample_data = sample_data.sort_values('eps_ax').reset_index(drop=True)
        
        # Get orientation
        orientation = sample_data['orientation'].iloc[0]
        
        # Determine analysis range (up to max_stress_fraction of max stress)
        max_stress = sample_data['sigma_d_psi'].max()
        stress_threshold = max_stress * self.max_stress_fraction
        analysis_mask = sample_data['sigma_d_psi'] <= stress_threshold
        
        # Extract relevant columns
        eps_ax = sample_data['eps_ax'].values
        sigma_d = sample_data['sigma_d_psi'].values
        eps_rad_a = sample_data['eps_rad_a'].values
        eps_rad_b = sample_data['eps_rad_b'].values
        
        # Smooth the data
        eps_ax_smooth = self._smooth_data(eps_ax)
        sigma_d_smooth = self._smooth_data(sigma_d)
        eps_rad_a_smooth = self._smooth_data(eps_rad_a)
        eps_rad_b_smooth = self._smooth_data(eps_rad_b)
        
        # Compute Young's Modulus (derivative of stress vs strain)
        YM = self._compute_derivative(sigma_d_smooth, eps_ax_smooth)
        
        # Compute Poisson's Ratios (negative derivative of radial vs axial strain)
        PRA = -self._compute_derivative(eps_rad_a_smooth, eps_ax_smooth)
        PRB = -self._compute_derivative(eps_rad_b_smooth, eps_ax_smooth)
        
        # Average radial strain for vertical samples
        if orientation == 'Vertical':
            eps_rad_avg = (eps_rad_a + eps_rad_b) / 2
            eps_rad_avg_smooth = self._smooth_data(eps_rad_avg)
            PRAvg = -self._compute_derivative(eps_rad_avg_smooth, eps_ax_smooth)
        else:
            PRAvg = np.full_like(PRA, np.nan)
        
        # Store processed data
        results = {
            'sample_id': sample_id,
            'orientation': orientation,
            'max_stress': max_stress,
            'stress_threshold': stress_threshold,
            'n_points_total': len(sample_data),
            'n_points_analysis': analysis_mask.sum(),
            
            # Full profiles
            'eps_ax': eps_ax,
            'sigma_d': sigma_d,
            'eps_rad_a': eps_rad_a,
            'eps_rad_b': eps_rad_b,
            'YM': YM,
            'PRA': PRA,
            'PRB': PRB,
            'PRAvg': PRAvg,
            'analysis_mask': analysis_mask,
            
            # Statistics within analysis range
            'YM_mean': np.nanmean(YM[analysis_mask]),
            'YM_std': np.nanstd(YM[analysis_mask]),
            'YM_median': np.nanmedian(YM[analysis_mask]),
            'YM_min': np.nanmin(YM[analysis_mask]),
            'YM_max': np.nanmax(YM[analysis_mask]),
            
            'PRA_mean': np.nanmean(PRA[analysis_mask]),
            'PRA_std': np.nanstd(PRA[analysis_mask]),
            'PRA_median': np.nanmedian(PRA[analysis_mask]),
            
            'PRB_mean': np.nanmean(PRB[analysis_mask]),
            'PRB_std': np.nanstd(PRB[analysis_mask]),
            'PRB_median': np.nanmedian(PRB[analysis_mask]),
        }
        
        if orientation == 'Vertical':
            results['PRAvg_mean'] = np.nanmean(PRAvg[analysis_mask])
            results['PRAvg_std'] = np.nanstd(PRAvg[analysis_mask])
            results['PRAvg_median'] = np.nanmedian(PRAvg[analysis_mask])
        
        self.results[sample_id] = results
        return results
    
    def process_all_samples(self):
        """Process all samples in the dataset"""
        sample_ids = self.df['sample_id'].unique()
        print(f"\nProcessing {len(sample_ids)} samples...")
        
        for i, sample_id in enumerate(sample_ids, 1):
            print(f"  [{i}/{len(sample_ids)}] {sample_id}")
            self.process_sample(sample_id)
        
        print("\nProcessing complete!")
        return self.results
    
    def get_summary_statistics(self):
        """Generate summary statistics table for all samples"""
        if not self.results:
            self.process_all_samples()
        
        summary_data = []
        for sample_id, res in self.results.items():
            row = {
                'Sample ID': sample_id,
                'Orientation': res['orientation'],
                'Max Stress (psi)': res['max_stress'],
                'Analysis Points': res['n_points_analysis'],
                'YM Mean (psi)': res['YM_mean'],
                'YM Std (psi)': res['YM_std'],
                'YM Median (psi)': res['YM_median'],
                'PRA Mean': res['PRA_mean'],
                'PRA Std': res['PRA_std'],
                'PRB Mean': res['PRB_mean'],
                'PRB Std': res['PRB_std'],
            }
            
            if res['orientation'] == 'Vertical':
                row['PRAvg Mean'] = res.get('PRAvg_mean', np.nan)
                row['PRAvg Std'] = res.get('PRAvg_std', np.nan)
            
            summary_data.append(row)
        
        return pd.DataFrame(summary_data)
    
    def plot_sample(self, sample_id, save_path=None):
        """
        Create comprehensive plots for a single sample
        
        Parameters:
        -----------
        sample_id : str
            Sample identifier
        save_path : str, optional
            Path to save figure
        """
        if sample_id not in self.results:
            self.process_sample(sample_id)
        
        res = self.results[sample_id]
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle(f'Sample: {sample_id} ({res["orientation"]})', fontsize=14, fontweight='bold')
        
        # Stress-strain curve
        ax = axes[0, 0]
        ax.plot(res['eps_ax'], res['sigma_d'], 'b-', linewidth=1.5, label='Full curve')
        ax.axvline(res['eps_ax'][res['analysis_mask']][-1], color='r', linestyle='--', 
                   label=f'Analysis limit ({self.max_stress_fraction*100:.0f}% max stress)')
        ax.set_xlabel('Axial Strain', fontsize=11)
        ax.set_ylabel('Deviatoric Stress (psi)', fontsize=11)
        ax.set_title('Stress-Strain Curve', fontsize=12)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Young's Modulus profile
        ax = axes[0, 1]
        mask = res['analysis_mask']
        ax.plot(res['eps_ax'][mask], res['YM'][mask], 'g-', linewidth=1.5)
        ax.axhline(res['YM_mean'], color='r', linestyle='--', label=f'Mean: {res["YM_mean"]:.0f} psi')
        ax.fill_between(res['eps_ax'][mask], 
                        res['YM_mean'] - res['YM_std'], 
                        res['YM_mean'] + res['YM_std'],
                        alpha=0.2, color='r', label=f'±1 Std: {res["YM_std"]:.0f} psi')
        ax.set_xlabel('Axial Strain', fontsize=11)
        ax.set_ylabel("Young's Modulus (psi)", fontsize=11)
        ax.set_title("Young's Modulus Evolution", fontsize=12)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Poisson's Ratio A profile
        ax = axes[1, 0]
        ax.plot(res['eps_ax'][mask], res['PRA'][mask], 'purple', linewidth=1.5, label='PRA')
        ax.axhline(res['PRA_mean'], color='r', linestyle='--', label=f'Mean: {res["PRA_mean"]:.3f}')
        ax.set_xlabel('Axial Strain', fontsize=11)
        ax.set_ylabel("Poisson's Ratio A", fontsize=11)
        ax.set_title("Poisson's Ratio A Evolution", fontsize=12)
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_ylim([0, 1])
        
        # Poisson's Ratio B (or Avg for vertical)
        ax = axes[1, 1]
        if res['orientation'] == 'Vertical':
            ax.plot(res['eps_ax'][mask], res['PRAvg'][mask], 'orange', linewidth=1.5, label='PRAvg')
            ax.axhline(res['PRAvg_mean'], color='r', linestyle='--', 
                      label=f'Mean: {res["PRAvg_mean"]:.3f}')
            ax.set_title("Poisson's Ratio Average Evolution", fontsize=12)
        else:
            ax.plot(res['eps_ax'][mask], res['PRB'][mask], 'brown', linewidth=1.5, label='PRB')
            ax.axhline(res['PRB_mean'], color='r', linestyle='--', 
                      label=f'Mean: {res["PRB_mean"]:.3f}')
            ax.set_title("Poisson's Ratio B Evolution", fontsize=12)
        
        ax.set_xlabel('Axial Strain', fontsize=11)
        ax.set_ylabel("Poisson's Ratio", fontsize=11)
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_ylim([0, 1])
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Saved plot to {save_path}")
        
        plt.show()
    
    def plot_comparison(self, sample_ids=None, save_path=None):
        """
        Create comparison plots across multiple samples
        
        Parameters:
        -----------
        sample_ids : list, optional
            List of sample IDs to compare. If None, uses all samples
        save_path : str, optional
            Path to save figure
        """
        if not self.results:
            self.process_all_samples()
        
        if sample_ids is None:
            sample_ids = list(self.results.keys())
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('Multi-Sample Comparison', fontsize=14, fontweight='bold')
        
        colors = plt.cm.tab10(np.linspace(0, 1, len(sample_ids)))
        
        for i, sample_id in enumerate(sample_ids):
            res = self.results[sample_id]
            mask = res['analysis_mask']
            color = colors[i]
            
            # Stress-strain curves
            axes[0, 0].plot(res['eps_ax'], res['sigma_d'], color=color, 
                           linewidth=1, alpha=0.7, label=sample_id)
            
            # YM profiles
            axes[0, 1].plot(res['eps_ax'][mask], res['YM'][mask], color=color,
                           linewidth=1, alpha=0.7)
            
            # PRA profiles
            axes[1, 0].plot(res['eps_ax'][mask], res['PRA'][mask], color=color,
                           linewidth=1, alpha=0.7)
            
            # PRB or PRAvg profiles
            if res['orientation'] == 'Vertical':
                axes[1, 1].plot(res['eps_ax'][mask], res['PRAvg'][mask], color=color,
                               linewidth=1, alpha=0.7, label=sample_id)
            else:
                axes[1, 1].plot(res['eps_ax'][mask], res['PRB'][mask], color=color,
                               linewidth=1, alpha=0.7, label=sample_id)
        
        axes[0, 0].set_xlabel('Axial Strain')
        axes[0, 0].set_ylabel('Deviatoric Stress (psi)')
        axes[0, 0].set_title('Stress-Strain Curves')
        axes[0, 0].legend(fontsize=8, loc='best')
        axes[0, 0].grid(True, alpha=0.3)
        
        axes[0, 1].set_xlabel('Axial Strain')
        axes[0, 1].set_ylabel("Young's Modulus (psi)")
        axes[0, 1].set_title("Young's Modulus Profiles")
        axes[0, 1].grid(True, alpha=0.3)
        
        axes[1, 0].set_xlabel('Axial Strain')
        axes[1, 0].set_ylabel("Poisson's Ratio A")
        axes[1, 0].set_title("Poisson's Ratio A Profiles")
        axes[1, 0].grid(True, alpha=0.3)
        axes[1, 0].set_ylim([0, 1])
        
        axes[1, 1].set_xlabel('Axial Strain')
        axes[1, 1].set_ylabel("Poisson's Ratio")
        axes[1, 1].set_title("Poisson's Ratio B/Avg Profiles")
        axes[1, 1].grid(True, alpha=0.3)
        axes[1, 1].set_ylim([0, 1])
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Saved comparison plot to {save_path}")
        
        plt.show()
    
    def plot_statistics_summary(self, save_path=None):
        """
        Create box plots summarizing YM and PR statistics across all samples
        """
        if not self.results:
            self.process_all_samples()
        
        # Collect data
        ym_data = [res['YM_mean'] for res in self.results.values()]
        pra_data = [res['PRA_mean'] for res in self.results.values()]
        prb_data = [res['PRB_mean'] for res in self.results.values()]
        
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        fig.suptitle('Statistical Summary Across All Samples', fontsize=14, fontweight='bold')
        
        # YM box plot
        axes[0].boxplot(ym_data, labels=['All Samples'])
        axes[0].set_ylabel("Young's Modulus (psi)", fontsize=11)
        axes[0].set_title(f"YM Distribution (n={len(ym_data)})", fontsize=12)
        axes[0].grid(True, alpha=0.3, axis='y')
        
        # PRA box plot
        axes[1].boxplot(pra_data, labels=['All Samples'])
        axes[1].set_ylabel("Poisson's Ratio A", fontsize=11)
        axes[1].set_title(f"PRA Distribution (n={len(pra_data)})", fontsize=12)
        axes[1].grid(True, alpha=0.3, axis='y')
        
        # PRB box plot
        axes[2].boxplot(prb_data, labels=['All Samples'])
        axes[2].set_ylabel("Poisson's Ratio B", fontsize=11)
        axes[2].set_title(f"PRB Distribution (n={len(prb_data)})", fontsize=12)
        axes[2].grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Saved statistics plot to {save_path}")
        
        plt.show()


# Example usage
if __name__ == "__main__":
    # Initialize analyzer
    analyzer = TriaxialAnalyzer(
        filepath='your_data.csv',  # Replace with your file path
        max_stress_fraction=0.70,  # Analyze up to 70% of max stress
        smooth_window=11,          # Smoothing window size
        smooth_order=2             # Polynomial order for smoothing
    )
    
    # Process all samples
    analyzer.process_all_samples()
    
    # Get summary statistics
    summary = analyzer.get_summary_statistics()
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    print(summary.to_string(index=False))
    
    # Save summary to CSV
    summary.to_csv('triaxial_summary_statistics.csv', index=False)
    print("\nSaved summary statistics to 'triaxial_summary_statistics.csv'")
    
    # Plot individual samples (first 3 as examples)
    sample_ids = list(analyzer.results.keys())[:3]
    for sample_id in sample_ids:
        analyzer.plot_sample(sample_id, save_path=f'plot_{sample_id}.png')
    
    # Create comparison plot
    analyzer.plot_comparison(save_path='comparison_plot.png')
    
    # Create statistical summary
    analyzer.plot_statistics_summary(save_path='statistics_summary.png')