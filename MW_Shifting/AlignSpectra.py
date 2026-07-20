import numpy as np
import pandas as pd
from io import StringIO
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.optimize import curve_fit

try:
    import urQRd as denoise
except ImportError:
    denoise = None


class AlignSpectra:
    def __init__(self, file_path=None, is_old=True, a=None, b=None, columns=None, df=None,
                 shift_model=None, shift_params=None):
        """
        Initialize the processor.

        :param file_path: String - Path to spectra data file
        :param is_old: Boolean - Whether data is from old equipment (needs shifting)
        :param a: Integer/float - (Deprecated) Shift parameter for log model (a*ln(x)+b)
        :param b: Integer/float - (Deprecated) Shift parameter for log model (a*ln(x)+b)
        :param columns: Range - Spectral columns to use (default: 50-700)
        :param shift_model: String - Name of shift model ('log', 'linear', 'quadratic', 'exponential')
        :param shift_params: dict - Coefficients for the chosen model, e.g., {'a':1, 'b':2, 'c':3}
        """
        self.file_path = file_path
        self.is_old = is_old
        self.df = df
        self.header = ""
        self.processed_data = None
        self.columns = columns if columns else [str(a) for a in range(50, 701)]
        self.types = {col: int for col in self.columns}
        self.shift_params_calibrated = False

        # New shift model attributes
        self.shift_model_type = None
        self.shift_params = {}

        # Handle deprecated a,b arguments (log model)
        if a is not None or b is not None:
            if a is None or b is None:
                raise ValueError("Both 'a' and 'b' must be provided when using the log shift model.")
            warnings.warn("Using 'a' and 'b' is deprecated. Use shift_model='log' and shift_params={'a':a, 'b':b} instead.",
                          DeprecationWarning, stacklevel=2)
            self.shift_model_type = 'log'
            self.shift_params = {'a': a, 'b': b}
        else:
            # Use new explicit model if provided
            if shift_model is not None and shift_params is not None:
                self.shift_model_type = shift_model
                self.shift_params = shift_params

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

    def load_df(self, df):
        """
        Load in a dataframe

        :param df: DataFrame - Pandas DataFrame of the 3D file
        :return: Boolean - If the DataFrame loads successfully
        """
        try:
            self.df = df
            return True
        except:
            return False

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
        Apply shift function to scan numbers according to the chosen model.

        :param x: Integer/float - Original scan number
        :return: Shifted scan number
        """
        if not self.is_old or self.shift_model_type is None:
            return x
        if self.shift_params is None or not self.shift_params:
            raise ValueError("Shift parameters have not been set.")

        model = self.shift_model_type
        p = self.shift_params
        if model == 'log':
            if x <= 0:
                return x
            return x + (p['a'] * np.log(x) + p['b'])
        elif model == 'linear':
            return x + (p['a'] * x + p['b'])
        elif model == 'quadratic':
            return x + (p['a'] * x**2 + p['b'] * x + p['c'])
        elif model == 'exponential':
            return x + (p['a'] * np.exp(p['b'] * x) + p['c'])
        else:
            raise ValueError(f"Unknown shift model: '{model}'")

    def set_shift_parameters(self, model, params=None):
        """
        Set shift model and parameters.

        :param model: String or numeric - Model name ('log','linear','quadratic','exponential')
                     or, for backward compatibility, numeric 'a' when used as set_shift_parameters(a, b)
        :param params: dict or numeric - Model parameters as dict, or 'b' when model is numeric 'a'
        """
        # Backward compatibility: old signature set_shift_parameters(a, b)
        if isinstance(model, (int, float)) and params is not None and isinstance(params, (int, float)):
            warnings.warn("Using set_shift_parameters(a,b) is deprecated. Use set_shift_parameters('log', {'a':a,'b':b}).",
                          DeprecationWarning, stacklevel=2)
            self.shift_model_type = 'log'
            self.shift_params = {'a': model, 'b': params}
            self.shift_params_calibrated = True
            return

        # New usage: model is string, params is dict
        if not isinstance(model, str) or not isinstance(params, dict):
            raise ValueError("set_shift_parameters expects a model name (str) and a params dictionary.")

        valid_models = ['log', 'linear', 'quadratic', 'exponential']
        if model not in valid_models:
            raise ValueError(f"Model must be one of {valid_models}")
        # Check required keys
        required = {'log': ['a', 'b'], 'linear': ['a', 'b'],
                    'quadratic': ['a', 'b', 'c'], 'exponential': ['a', 'b', 'c']}
        for key in required[model]:
            if key not in params:
                raise ValueError(f"Parameter '{key}' missing for model '{model}'")

        self.shift_model_type = model
        self.shift_params = params
        self.shift_params_calibrated = True

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
                                   threshold=5, peak_threshold=18000, verbose=False,
                                   models=None):
        """
        Calibrate optimal shift parameters across multiple systematic shift models.

        :param file_pairs: List of (old_file_path, new_file_path) tuples
        :param window_step: Window size increment for optimization
        :param max_window: Maximum window size to test
        :param threshold: Error threshold for accuracy calculation
        :param peak_threshold: Minimum intensity for peak detection
        :param verbose: Whether to print progress
        :param models: List of model names to try (default: ['log','linear','quadratic','exponential'])
        :return: (best_model_name, best_params_dict, best_window, best_accuracy)
        """
        if models is None:
            models = ['log', 'linear', 'quadratic', 'exponential']

        best_global_accuracy = -1
        best_global_window = None
        best_global_model = None
        best_global_params = None

        for win in range(window_step, max_window + 1, window_step):
            # Collect all clean (x, offset) pairs across all file pairs for this window
            all_x = []
            all_offset = []

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
                    lower, upper = q1 - 1.5 * iqr, q3 + 1.5 * iqr
                    mask = (y_arr >= lower) & (y_arr <= upper)
                    clean_x, clean_y = x_arr[mask], y_arr[mask]

                    if len(clean_x) < 3:
                        continue

                    all_x.extend(clean_x.tolist())
                    all_offset.extend(clean_y.tolist())

                except Exception as e:
                    if verbose:
                        print(f"Error processing {old_path} and {new_path}: {str(e)}")
                    continue

            if len(all_x) < 3:
                if verbose:
                    print(f"Window {win}: insufficient data")
                continue

            x_arr = np.array(all_x)
            y_arr = np.array(all_offset)

            # Try each model on the aggregated data for this window
            for model in models:
                try:
                    if model == 'log':
                        # offset = a*ln(x) + b
                        log_x = np.log(x_arr)
                        coeffs = np.polyfit(log_x, y_arr, 1)  # [a, b]
                        a, b = coeffs[0], coeffs[1]
                        predicted = a * log_x + b
                        params = {'a': a, 'b': b}
                    elif model == 'linear':
                        coeffs = np.polyfit(x_arr, y_arr, 1)
                        a, b = coeffs[0], coeffs[1]
                        predicted = a * x_arr + b
                        params = {'a': a, 'b': b}
                    elif model == 'quadratic':
                        coeffs = np.polyfit(x_arr, y_arr, 2)
                        a, b, c = coeffs[0], coeffs[1], coeffs[2]
                        predicted = a * x_arr**2 + b * x_arr + c
                        params = {'a': a, 'b': b, 'c': c}
                    elif model == 'exponential':
                        # offset = a*exp(b*x) + c
                        def exp_func(x, a, b, c):
                            return a * np.exp(b * x) + c
                        # Initial guess: small negative exponent, magnitude around median offset
                        init_a = np.median(y_arr) * 0.5
                        init_b = -1e-4
                        init_c = np.median(y_arr) * 0.5
                        popt, _ = curve_fit(exp_func, x_arr, y_arr,
                                            p0=[init_a, init_b, init_c],
                                            maxfev=5000)
                        predicted = exp_func(x_arr, *popt)
                        params = {'a': popt[0], 'b': popt[1], 'c': popt[2]}
                    else:
                        continue

                    # Accuracy: fraction of predictions within 'threshold' of actual offset
                    errors = np.abs(predicted - y_arr)
                    good = np.sum(errors <= threshold)
                    accuracy = good / len(y_arr)

                    if verbose:
                        print(f"Window {win}, Model {model}: accuracy={accuracy:.4f}, params={params}")

                    if accuracy > best_global_accuracy:
                        best_global_accuracy = accuracy
                        best_global_window = win
                        best_global_model = model
                        best_global_params = params

                except Exception as e:
                    if verbose:
                        print(f"Window {win}, Model {model}: fitting failed ({e})")
                    continue

        if best_global_model is None:
            raise RuntimeError("Calibration failed - no valid model found for any window")

        if verbose:
            print(f"\nOptimal model: {best_global_model}")
            print(f"Optimal window: {best_global_window}")
            print(f"Optimal parameters: {best_global_params}")
            print(f"Accuracy: {best_global_accuracy*100:.2f}%")

        return best_global_model, best_global_params, best_global_window, best_global_accuracy
