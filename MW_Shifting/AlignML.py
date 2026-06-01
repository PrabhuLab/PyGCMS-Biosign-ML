import numpy as np
import pandas as pd
from io import StringIO
import os
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
try:
    import urQRd as denoise
except ImportError:
    denoise = None

class AlignSpectra:
    def __init__(self, file_path=None, is_old=True, a=None, b=None, columns=None):
        """
        Initialize the SpectraProcessor.
        
        :param file_path: String - Path to spectra data file
        :param is_old: Boolean - Whether data is from old equipment (needs shifting)
        :param a: Integer - Shift parameter a (log coefficient) (a * ln(x) + b)
        :param b: Integer - Shift parameter b (constant offset) (a * ln(x) + b)
        :param columns: Range - Spectral columns to use (default: 50-700)
        """
        self.file_path = file_path
        self.is_old = is_old
        self.a = a
        self.b = b
        self.df = None
        self.header = ""
        self.processed_data = None
        self.columns = columns if columns else [str(a) for a in range(50, 701)]
        self.types = {col: int for col in self.columns}
        self.shift_params_calibrated = False
        
        if file_path:
            self.load_file()
    
    def load_file(self, file_path=None):
        """
        Load and parse spectra data file.
        
        :param file_path: String - Optional override of file path
        :return: Loaded DataFrame
        """
        path = file_path or self.file_path
        if not path:
            raise ValueError("No file path provided")
            
        with open(path, 'r') as f:
            content = f.read()
        
        if "Scan," not in content:
            raise ValueError("File format error: 'Scan,' header not found")
        
        parts = content.split("Scan,", 1)
        self.header = parts[0]
        data_part = "Scan," + parts[1].replace(" ", "")
        
        self.df = pd.read_csv(
            StringIO(data_part), 
            dtype=self.types, 
            index_col=[0]
        ).loc[:, lambda df: ~df.columns.str.contains('^Unnamed')]
        
        return self.df
    
    @staticmethod
    def baseline(y, lam=1e5, p=0.01, n_iter=10):
        """
        Calculate baseline using asymmetric least squares.
        
        :param y: List - Raw intensity values
        :param lam: Integer - Smoothness parameter
        :param p: Float - Asymmetry parameter
        :param n_iter: Integer - Number of iterations
        :return: Baseline estimate
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
    
    def preprocess(self, df=None, columns=None):
        """
        Preprocess spectra data: baseline correction and denoising.
        
        :param df: Dataframe - Optional DataFrame to process
        :param columns: List - Spectral columns to use
        :return: Processed data array
        """
        if denoise is None:
            raise ImportError("urQRd package is required for preprocessing")
            
        df = df if df is not None else self.df
        columns = columns or self.columns
        
        if df is None:
            raise ValueError("No data available for preprocessing")
            
        intensity = np.sum(df[columns], axis=1).to_numpy()
        base = self.baseline(intensity)
        corrected = intensity - base
        denoised = denoise.urQRd(
            corrected.astype(complex), 
            k=100, 
            orda=min(400, len(corrected) // 2), 
            iterations=1
        )
        self.processed_data = np.real(denoised)
        return self.processed_data
    
    @staticmethod
    def detect_maxima(y_vals, x_vals, step=110, threshold=18000):
        """
        Identify peak maxima in processed data.
        
        :param y_vals: List - Processed intensity values
        :param x_vals: List - Scan numbers
        :param step: Integer - Window size for peak detection
        :param threshold: Integer - Minimum intensity for peaks
        :return: List of peak positions
        """
        peaks_x = []
        for i in range(0, len(y_vals), step):
            window_y = y_vals[i:i + step]
            window_x = x_vals[i:i + step]
            if len(window_y) == 0:
                continue
            idxs = np.argsort(window_y)[-2:][::-1]
            for idx in idxs:
                if window_y[idx] > threshold:
                    peaks_x.append(window_x[idx])
        return peaks_x
    
    def shift_function(self, x):
        """
        Apply shift function to scan numbers.
        
        :param x: Integer - Original scan number
        :return: Shifted scan number
        """
        if self.a == None:
            raise ValueError("No valid \"a\" value. Maybe try calibrating the shift parameters (calibrate_shift_parameters) first and passing the \"a\" and \"b\" in the SpectraProcessor() class")
        if self.b == None:
            raise ValueError("No valid \"b\" value. Maybe try calibrating the shift parameters (calibrate_shift_parameters) first and passing the \"a\" and \"b\" in the SpectraProcessor() class")
        if x <= 0:
            return x
        return x + (self.a * np.log(x) + self.b)
    
    def save_shifted(self, output_path):
        """
        Save shifted data to a new file.
        
        :param output_path: String - Output file path
        :return: Output file path
        """
        if self.df is None:
            raise ValueError("No data available for shifting")
            
        if self.is_old:
            shifted_index = self.df.index.map(self.shift_function)
            self.df.index = shifted_index
            
        df_out = self.df.reset_index().rename(columns={'index': 'Scan'})
        new_csv = df_out.to_csv(index=False)
        new_content = self.header + new_csv
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w') as f:
            f.write(new_content)
            
        return output_path
    
    def plot_spectra(self, processed=True, title=None):
        """
        Plot spectra data.
        
        :param processed: Boolean - Whether to plot processed data
        :param title: String - Plot title
        :return: Matplotlib figure
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
            y = np.sum(self.df[self.columns], axis=1)
            label = "Raw Data"
            
        ax.plot(self.df.index, y, label=label)
        ax.set_xlabel("Scan Number")
        ax.set_ylabel("Intensity")
        title = title or f"Spectra: {os.path.basename(self.file_path)}"
        ax.set_title(title)
        ax.legend()
        ax.grid(True)
        
        return fig
    
    @classmethod
    def calibrate_shift_parameters(cls, file_pairs, window_step=5, max_window=200, 
                                   threshold=5, peak_threshold=18000, verbose=False):
        """
        Calibrate optimal shift parameters (a, b) by comparing old and new equipment files.
        
        :param file_pairs: List - List of (old_file_path, new_file_path) tuples
        :param window_step: Integer - Window size increment for optimization
        :param max_window: Integer - Maximum window size to test
        :param threshold: Integer - Error threshold for accuracy calculation
        :param peak_threshold: Integer - Minimum intensity for peak detection
        :param verbose: Boolean - Whether to print progress
        :return: (a, b, best_window, best_accuracy)
        """
        best_window = None
        best_accuracy = 0
        best_model = (0, 0)
        all_logx = []
        all_y = []
        
        for win in range(window_step, max_window + 1, window_step):
            total_good = 0
            total_count = 0
            win_logx = []
            win_y = []
            
            for old_path, new_path in file_pairs:
                try:
                    proc_old = cls(old_path)
                    df_old = proc_old.df
                    processed_old = proc_old.preprocess()
                    x_old = df_old.index.tolist()
                    peaks_old = cls.detect_maxima(processed_old, x_old, win, peak_threshold)
                    
                    proc_new = cls(new_path)
                    df_new = proc_new.df
                    processed_new = proc_new.preprocess()
                    x_new = df_new.index.tolist()
                    peaks_new = cls.detect_maxima(processed_new, x_new, win, peak_threshold)
                    
                    n = min(len(peaks_old), len(peaks_new))
                    if n < 3:
                        continue
                    
                    offsets = []
                    peak_positions = []
                    for i in range(n):
                        nx = peaks_new[i]
                        ox = peaks_old[i]
                        offset = ox - nx
                        if 0 <= offset <= 35 and nx <= 4000:
                            offsets.append(offset)
                            peak_positions.append(nx)
                    
                    if len(offsets) < 3:
                        continue
                    
                    x_arr, y_arr = np.array(peak_positions), np.array(offsets)
                    q1, q3 = np.percentile(y_arr, [25, 75])
                    iqr = q3 - q1
                    lower, upper = q1 - 1.5*iqr, q3 + 1.5*iqr
                    mask = (y_arr >= lower) & (y_arr <= upper)
                    clean_x, clean_y = x_arr[mask], y_arr[mask]
                    
                    if len(clean_x) < 3:
                        continue
                    
                    log_x = np.log(clean_x)
                    win_logx.extend(log_x)
                    win_y.extend(clean_y)
                    
                except Exception as e:
                    if verbose:
                        print(f"Error processing {old_path} and {new_path}: {str(e)}")
                    continue
            
            if len(win_logx) < 3:
                if verbose:
                    print(f"Window {win}: insufficient data")
                continue
            
            a_win, b_win = np.polyfit(win_logx, win_y, 1)
            
            predicted = np.array(win_logx) * a_win + b_win
            errors = np.abs(predicted - np.array(win_y))
            good = np.sum(errors <= threshold)
            count = len(errors)
            accuracy = good / count if count else 0
            
            if verbose:
                print(f"Window {win}: accuracy={accuracy:.4f}, a={a_win:.4f}, b={b_win:.4f}")
            
            if accuracy > best_accuracy:
                best_accuracy = accuracy
                best_window = win
                best_model = (a_win, b_win)
                all_logx = win_logx
                all_y = win_y
        
        if best_model == (0, 0):
            raise RuntimeError("Calibration failed - no valid models found")
        
        a, b = best_model
        if verbose:
            print(f"Optimal window: {best_window}")
            print(f"Optimal model: {a:.5f}*ln(x) + {b:.5f}")
            print(f"Accuracy: {best_accuracy*100:.2f}%")
        
        return a, b, best_window, best_accuracy
    
    def set_shift_parameters(self, a, b):
        """
        Set shift parameters after calibration. Function: a * ln(x) + b
        
        :param a: Float - The "a" value in the function
        :param b: Float - The "b" value in the function
        """
        self.a = a
        self.b = b
        self.shift_params_calibrated = True
