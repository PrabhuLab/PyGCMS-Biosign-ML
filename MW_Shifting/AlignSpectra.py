import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.interpolate import interp1d

try:
    import urQRd as denoise
except ImportError:
    denoise = None

class AlignSpectra:
    def __init__(self, file_path=None, df=None, is_old=True, shift_type='logarithmic', 
                 a=None, b=None, c=None, y_cols=None, x_label="X-axis", y_label="Intensity"):
        """
        Initialize the SpectraProcessor.
        
        :param file_path: String - Path to spectra data file (CSV)
        :param df: DataFrame - Optional pre-loaded pandas DataFrame
        :param is_old: Boolean - Whether data needs shifting
        :param shift_type: String - Type of shift ('linear', 'logarithmic', 'exponential', 'polynomial', 'power')
        :param a, b, c: Floats - Shift parameters
        :param y_cols: List - Specific columns to sum (for 2D data). If None, sums all numeric columns.
        :param x_label: String - Label for the X-axis (e.g., 'Wavenumber', 'm/z', 'Retention Time')
        :param y_label: String - Label for the Y-axis (e.g., 'Absorbance', 'Intensity')
        """
        self.file_path = file_path
        self.df = df
        self.is_old = is_old
        self.shift_type = shift_type
        self.a = a
        self.b = b
        self.c = c
        self.y_cols = y_cols
        self.x_label = x_label
        self.y_label = y_label
        self.processed_data = None
        self.shift_params_calibrated = False
        
        if self.df is None and self.file_path:
            self.load_file()
    
    def load_file(self, file_path=None, index_col=0, skiprows=0):
        """
        Load spectra data from a generic CSV.
        
        :param file_path: String - Optional override of file path
        :param index_col: Integer - Column to use as the X-axis
        :param skiprows: Integer - Number of metadata rows to skip at the top
        :return: Loaded DataFrame
        """
        path = file_path or self.file_path
        if not path:
            raise ValueError("No file path provided")
            
        self.df = pd.read_csv(path, index_col=index_col, skiprows=skiprows)
        # Drop empty unnamed columns often left by trailing commas
        self.df = self.df.loc[:, ~self.df.columns.str.contains('^Unnamed')]
        
        return self.df
    
    @staticmethod
    def baseline(y, lam=1e5, p=0.01, n_iter=10):
        """
        Calculate baseline using asymmetric least squares.
        Universally applicable to MS, Raman, IR, XRD, etc.
        """
        m = len(y)
        diagonals = [np.ones(m), -2 * np.ones(m), np.ones(m)]
        D = diags(diagonals, [0, 1, 2], shape=(m - 2, m), format='csr')
        DTD = D.T @ D
        w = np.ones(m)
        for _ in range(n_iter):
            W = diags(w, 0)
            Z = W + lam * DTD
            z = spsolve(Z, w * y)
            w = p * (y > z) + (1 - p) * (y < z)
        return z
    
    def preprocess(self, df=None, y_cols=None):
        """
        Preprocess spectra data: baseline correction and optional denoising.
        """
        df = df if df is not None else self.df
        y_cols = y_cols or self.y_cols
        
        if df is None:
            raise ValueError("No data available for preprocessing")
            
        # Handle 1D or 2D data
        if y_cols is not None:
            intensity = np.sum(df[y_cols], axis=1).to_numpy()
        else:
            intensity = df.sum(axis=1).to_numpy()
            
        base = self.baseline(intensity)
        corrected = intensity - base
        
        if denoise is not None:
            denoised = denoise.urQRd(
                corrected.astype(complex), 
                k=100, 
                orda=min(400, max(1, len(corrected) // 2)), 
                iterations=1
            )
            self.processed_data = np.real(denoised)
        else:
            self.processed_data = corrected
            
        return self.processed_data
    
    @staticmethod
    def detect_maxima(y_vals, x_vals, prominence=1000, distance=10):
        """
        Identify peak maxima using SciPy's generalized peak finder.
        
        :param y_vals: Array - Processed intensity values
        :param x_vals: Array - X-axis values
        :param prominence: Float - Minimum height of peak relative to baseline
        :param distance: Integer - Minimum horizontal distance between peaks
        :return: Array of peak X-positions
        """
        peaks, _ = find_peaks(y_vals, prominence=prominence, distance=distance)
        return np.array(x_vals)[peaks]
    
    def shift_function(self, x):
        """
        Apply shift function to X-axis values.
        """
        if self.a is None or self.b is None:
            raise ValueError("Shift parameters missing. Set them manually or run calibrate_shift_parameters().")
        if self.shift_type in ['exponential', 'polynomial', 'power'] and self.c is None:
            raise ValueError(f"The '{self.shift_type}' shift requires a 'c' parameter.")
            
        if self.shift_type == 'logarithmic':
            # Protect against log(0) or negative logs
            return x if x <= 0 else x + (self.a * np.log(x) + self.b)
        elif self.shift_type == 'linear':
            return x + (self.a * x + self.b)
        elif self.shift_type == 'exponential':
            return x + (self.a * np.exp(self.b * x) + self.c)
        elif self.shift_type == 'polynomial':
            return x + (self.a * (x ** 2) + self.b * x + self.c)
        elif self.shift_type == 'power':
            return x if x <= 0 else x + (self.a * (x ** self.b) + self.c)
        else:
            raise ValueError(f"Unknown shift type: {self.shift_type}")
    
    def apply_shift(self, interpolate=True):
        """
        Applies mathematical shift to the X-axis and optionally interpolates 
        data back onto the original X-grid for proper alignment.
        """
        if self.df is None:
            raise ValueError("No data available for shifting")
            
        original_x = self.df.index.to_numpy(dtype=float)
        shifted_x = np.array([self.shift_function(x) for x in original_x])
        
        if interpolate:
            new_df = pd.DataFrame(index=original_x, columns=self.df.columns)
            for col in self.df.columns:
                # Interpolate Y onto original X grid. Data outside bounds becomes 0.
                f = interp1d(shifted_x, self.df[col], kind='linear', bounds_error=False, fill_value=0)
                new_df[col] = f(original_x)
            self.df = new_df
        else:
            self.df.index = shifted_x
            
        return self.df
    
    def save_shifted(self, output_path, interpolate=True):
        """
        Apply shifts and save data to a new file.
        """
        if self.is_old:
            self.apply_shift(interpolate=interpolate)
            
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        self.df.to_csv(output_path, index_label=self.x_label)
        return output_path
    
    def plot_spectra(self, processed=True, title=None):
        """
        Plot spectra data dynamically.
        """
        if self.df is None:
            raise ValueError("No data available for plotting")
            
        fig, ax = plt.subplots(figsize=(10, 6))
        
        if processed:
            if self.processed_data is None:
                self.preprocess()
            y = self.processed_data
            label = "Processed Data"
        else:
            y = self.df.sum(axis=1) if self.y_cols is None else self.df[self.y_cols].sum(axis=1)
            label = "Raw Data"
            
        ax.plot(self.df.index, y, label=label)
        ax.set_xlabel(self.x_label)
        ax.set_ylabel(self.y_label)
        title = title or ("Spectra Plot" if not self.file_path else f"Spectra: {os.path.basename(self.file_path)}")
        ax.set_title(title)
        ax.legend()
        ax.grid(True)
        
        return fig
    
    @classmethod
    def calibrate_shift_parameters(cls, file_pairs, dist_step=5, max_dist=100, 
                                   prominence=1000, error_threshold=5, max_offset=50, 
                                   shift_type='all', verbose=False):
        """
        Calibrate optimal shift parameters by comparing old and new equipment files.
        Uses closest-peak matching.
        
        :param file_pairs: List - List of (old_file_path, new_file_path) tuples
        :param dist_step: Integer - Peak distance sweep step
        :param max_dist: Integer - Maximum peak distance for sweep
        :param prominence: Float - Minimum peak height
        :param error_threshold: Float - Max allowable error for accuracy scoring
        :param max_offset: Float - Maximum allowable physical shift between matching peaks
        :param shift_type: String - Specific shift type to force, or 'all'
        :param verbose: Boolean - Whether to print progress
        :return: (best_type, best_model, best_dist, best_accuracy)
        """
        best_dist = None
        best_accuracy = 0
        best_model = None
        best_type = None
        
        types_to_test = ['linear', 'logarithmic', 'exponential', 'polynomial', 'power'] if shift_type == 'all' else [shift_type]
        
        for dist in range(dist_step, max_dist + 1, dist_step):
            win_x = []
            win_y = []
            
            for old_path, new_path in file_pairs:
                try:
                    proc_old = cls(old_path)
                    processed_old = proc_old.preprocess()
                    x_old = proc_old.df.index.to_numpy(dtype=float)
                    peaks_old = cls.detect_maxima(processed_old, x_old, prominence, dist)
                    
                    proc_new = cls(new_path)
                    processed_new = proc_new.preprocess()
                    x_new = proc_new.df.index.to_numpy(dtype=float)
                    peaks_new = cls.detect_maxima(processed_new, x_new, prominence, dist)
                    
                    if len(peaks_old) == 0 or len(peaks_new) == 0:
                        continue
                        
                    # Closest peak matching
                    for nx in peaks_new:
                        diffs = np.abs(peaks_old - nx)
                        min_idx = np.argmin(diffs)
                        offset = peaks_old[min_idx] - nx
                        
                        # Only keep matches that fall within realistic limits
                        if np.abs(offset) <= max_offset:
                            win_x.append(nx)
                            win_y.append(offset)
                            
                except Exception as e:
                    if verbose:
                        print(f"Error processing {old_path} and {new_path}: {str(e)}")
                    continue
            
            if len(win_x) < 3:
                continue
                
            x_arr = np.array(win_x)
            y_arr = np.array(win_y)
            
            # Filter outliers based on IQR
            q1, q3 = np.percentile(y_arr, [25, 75])
            iqr = q3 - q1
            mask = (y_arr >= q1 - 1.5*iqr) & (y_arr <= q3 + 1.5*iqr)
            x_arr, y_arr = x_arr[mask], y_arr[mask]
            
            if len(x_arr) < 3:
                continue
            
            for stype in types_to_test:
                try:
                    if stype == 'linear':
                        popt = np.polyfit(x_arr, y_arr, 1)
                        predicted = popt[0] * x_arr + popt[1]
                    elif stype == 'logarithmic':
                        # Ignore non-positive x values for log
                        pos_mask = x_arr > 0
                        popt = np.polyfit(np.log(x_arr[pos_mask]), y_arr[pos_mask], 1)
                        predicted = np.zeros_like(x_arr, dtype=float)
                        predicted[pos_mask] = popt[0] * np.log(x_arr[pos_mask]) + popt[1]
                    elif stype == 'polynomial':
                        popt = np.polyfit(x_arr, y_arr, 2)
                        predicted = popt[0] * (x_arr ** 2) + popt[1] * x_arr + popt[2]
                    elif stype == 'exponential':
                        def exp_func(x, a, b, c): return a * np.exp(b * x) + c
                        popt, _ = curve_fit(exp_func, x_arr, y_arr, maxfev=10000)
                        predicted = exp_func(x_arr, *popt)
                    elif stype == 'power':
                        def pow_func(x, a, b, c): return a * (x ** b) + c
                        pos_mask = x_arr > 0
                        popt, _ = curve_fit(pow_func, x_arr[pos_mask], y_arr[pos_mask], maxfev=10000)
                        predicted = np.zeros_like(x_arr, dtype=float)
                        predicted[pos_mask] = pow_func(x_arr[pos_mask], *popt)
                    else:
                        continue
                        
                    errors = np.abs(predicted - y_arr)
                    accuracy = np.sum(errors <= error_threshold) / len(errors)
                    
                    if verbose:
                        print(f"Dist {dist} [{stype}]: accuracy={accuracy:.4f}, params={popt}")
                    
                    if accuracy > best_accuracy:
                        best_accuracy = accuracy
                        best_dist = dist
                        best_model = tuple(popt)
                        best_type = stype
                        
                except Exception as e:
                    if verbose:
                        print(f"Fit failed for {stype} at distance {dist}: {e}")
                    continue
        
        if best_model is None:
            raise RuntimeError("Calibration failed - no valid models found. Try adjusting prominence or max_offset thresholds.")
        
        if verbose:
            print(f"\n--- Calibration Complete ---")
            print(f"Optimal peak distance: {best_dist}")
            print(f"Optimal shift type: {best_type}")
            print(f"Optimal model params: {best_model}")
            print(f"Accuracy: {best_accuracy*100:.2f}%")
        
        return best_type, best_model, best_dist, best_accuracy
    
    def set_shift_parameters(self, shift_type, a, b, c=None):
        """
        Set shift parameters after calibration.
        """
        self.shift_type = shift_type
        self.a = a
        self.b = b
        self.c = c
        self.shift_params_calibrated = True
