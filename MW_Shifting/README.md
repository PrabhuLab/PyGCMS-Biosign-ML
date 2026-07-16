# AlignSpectra

## Overview

Spectral and chromatographic data (such as py-GC-MS, Raman, FT-IR, and NMR) collected on different instruments often exhibit X-axis shifts, even when analyzing the same sample. These shifts can occur when transitioning between older and newer instruments, after routine maintenance, or due to laser/detector calibration differences, causing peaks to appear at different positions than before. 

This project aims to provide a universal framework for aligning datasets collected under different instrument conditions. The pipeline identifies corresponding peaks between reference datasets and models the relationship between their positions using best-fit mathematical functions (linear, logarithmic, exponential, polynomial, or power series). Once calibrated, the model can shift historical data and interpolate it back onto a standard grid so that peaks collected on different instruments align perfectly.

The processing pipeline includes baseline correction, optional denoising, peak detection, calibration of alignment parameters, interpolation, and export of corrected datasets. By aligning data collected before and after maintenance events or across multiple instruments, previously collected samples can be incorporated into the same analysis workflow without requiring manual adjustments.

---

## Features

* Import generic 1D and 2D spectral data from standard CSV files
* Handles discrete (Scan number) and continuous (Wavenumber, Retention Time) X-axes
* Baseline correction using Asymmetric Least Squares (ALS)
* Signal denoising using urQRd (optional)
* Automated peak detection using SciPy
* Calibration of spectral shifts between instruments
* Multi-model alignment (linear, logarithmic, exponential, polynomial, power)
* Data interpolation for standardized X-axis grids
* Visualization of raw and processed spectra
* Export of corrected datasets

---

## Installation

### Required Dependencies

```bash
pip install numpy pandas scipy matplotlib
```

### Optional Dependency

Preprocessing can optionally utilize urQRd for signal denoising, which should be located in this folder. Please note that this version of urQRd was modified by Miles Walters from the available copy of the original algorithm found [here](https://doi.org/10.1073/pnas.1306700111).   

Without urQRd, the pipeline will still execute baseline correction and alignment, but the specific denoising step will be bypassed.

---

## Data Format

The software operates on standard CSV datasets. It can process 1D single-trace spectra (e.g., a single Raman or FT-IR scan) or 2D matrices (e.g., GC-MS scan vs *m/z* ranges). 

Example CSV Structure:

```text
Instrument Metadata Header 1
Instrument Metadata Header 2
Scan,50,51,52,53,54,55,...
1,0,0,12,34,5,0,...
2,0,1,14,35,4,0,...
...
```

The loader automatically allows you to skip metadata rows and designate an X-axis column. 

---

## Basic Usage

### Loading a Dataset

```python
from AlignSpectra import AlignSpectra

processor = AlignSpectra(
    file_path="sample.csv",
    x_label="Scan",
    y_label="Intensity"
)
```

If your file has metadata headers to skip, you can load it manually:

```python
processor = AlignSpectra()

processor.load_file(
    file_path="sample.csv",
    index_col=0,
    skiprows=2
)
```

### Accessing the DataFrame

```python
df = processor.df

print(df.head())
```

---

## Preprocessing

Preprocessing consists of:

1. Summation of spectral intensities (if analyzing 2D matrix data)
2. Baseline correction
3. urQRd denoising (if installed)

```python
processed = processor.preprocess()
```

The processed signal is stored as a 1D NumPy array:

```python
processor.processed_data
```

---

## Plotting

### Processed Spectra

```python
import matplotlib.pyplot as plt

processor.preprocess()

processor.plot_spectra()

plt.show()
```

### Raw Spectra

```python
processor.plot_spectra(
    processed=False
)

plt.show()
```

---

## Peak Detection

Peak detection identifies local maxima in the processed spectra using SciPy's generalized peak finder.

```python
processed = processor.preprocess()

peaks = AlignSpectra.detect_maxima(
    processed,
    processor.df.index.tolist()
)

print(peaks)
```

### Custom Parameters

```python
peaks = AlignSpectra.detect_maxima(
    processed,
    processor.df.index.tolist(),
    prominence=1500,
    distance=20
)
```

| Parameter  | Description                                       |
| ---------- | ------------------------------------------------- |
| prominence | Minimum height of the peak relative to baseline   |
| distance   | Minimum horizontal X-axis distance between peaks  |

---

## Calibrating Alignment Parameters

Calibration requires paired datasets collected from different instruments or instrument states. The software can test all available shift models to find the mathematical function that aligns the data with the highest accuracy.

Example:

```python
pairs = [
    (
        "old_machine_run_1.csv",
        "new_machine_run_1.csv"
    ),
    (
        "old_machine_run_2.csv",
        "new_machine_run_2.csv"
    )
]

shift_type, model, dist, accuracy = (
    AlignSpectra.calibrate_shift_parameters(
        pairs,
        shift_type='all',
        verbose=True
    )
)

print("Optimal peak distance:", dist)
print("Optimal shift type:", shift_type)
print("Optimal model params:", model)
print(f"Accuracy: {accuracy * 100:.2f}%")
```

Output:

```text
Optimal peak distance: 15
Optimal shift type: polynomial
Optimal model params: (0.002, 1.45, 3.2)
Accuracy: 96.80%
```

### Calibration Method

For each paired dataset:

1. Preprocess both spectra
2. Detect peaks based on dynamic prominence
3. Match closest corresponding peaks within an allowable offset
4. Filter out mismatched outliers using interquartile ranges (IQR)
5. Fit models (linear, log, exponential, polynomial, power) to the peak shifts

The method iterates over peak `distance` parameters to determine the optimal density for creating the most accurate model. 

---

## Applying Alignment

After calibration, you can apply the model to an unaligned dataset:

```python
processor = AlignSpectra(
    file_path="old_sample.csv",
    is_old=True
)

processor.set_shift_parameters(
    shift_type=shift_type,
    a=model[0],
    b=model[1],
    c=model[2] if len(model) > 2 else None
)
```

The underlying data is mathematically shifted and, by default, interpolated back onto the original X-axis grid so datasets can be directly compared or subtracted.

---

## Exporting Aligned Data

```python
processor.save_shifted(
    "aligned/aligned_sample.csv",
    interpolate=True
)
```

The exported file:

* Applies mathematical X-axis correction
* Resamples intensities onto standard grid
* Writes a standard CSV file

---

## Complete Example

```python
from AlignSpectra import AlignSpectra

training_pairs = [
    (
        "old_reference.csv",
        "new_reference.csv"
    )
]

shift_type, model, _, accuracy = (
    AlignSpectra.calibrate_shift_parameters(
        training_pairs,
        shift_type='all',
        verbose=True
    )
)

print(
    f"Calibration accuracy: {accuracy:.2%}"
)

# Unpack parameters safely based on model size
a, b = model[0], model[1]
c = model[2] if len(model) > 2 else None

sample = AlignSpectra(
    file_path="old_unknown.csv",
    is_old=True
)

sample.set_shift_parameters(
    shift_type,
    a,
    b,
    c
)

sample.preprocess()

sample.plot_spectra()

sample.save_shifted(
    "aligned_unknown.csv"
)
```

---

## Constructor Parameters

```python
AlignSpectra(
    file_path=None,
    df=None,
    is_old=True,
    shift_type='logarithmic',
    a=None,
    b=None,
    c=None,
    y_cols=None,
    x_label="X-axis",
    y_label="Intensity"
)
```

| Parameter   | Description                                              |
| ----------- | -------------------------------------------------------- |
| file_path   | Path to spectral file (CSV)                              |
| df          | Optional pre-loaded pandas DataFrame                     |
| is_old      | Whether alignment should be applied                      |
| shift_type  | Mathematical shift ('linear', 'logarithmic', etc.)       |
| a, b, c     | Coefficients for the shift equation                      |
| y_cols      | Specific columns to sum (for 2D data). Blank sums all.   |
| x_label     | String label for the X-axis (e.g., 'Wavenumber')         |
| y_label     | String label for the Y-axis (e.g., 'Absorbance')         |

---

## Processing Workflow

```text
Raw Spectral Data --> Intensty Summation (Optional) --> Baseline Correction --> urQRd Denoising --> Peak Detection --> Calibration --> Model Selection --> Peak Alignment --> Grid Interpolation --> Corrected Dataset
```

---

## Error Handling

### Missing File

```python
ValueError:
No file path provided
```

### Missing Shift Parameters

```python
ValueError:
Shift parameters missing. Set them manually or run calibrate_shift_parameters().
```

### Missing Polynomial/Exponential Parameter

```python
ValueError:
The 'polynomial' shift requires a 'c' parameter.
```

### Calibration Failure

```python
RuntimeError:
Calibration failed - no valid models found. Try adjusting prominence or max_offset thresholds.
```

---

## Citation

If this software contributes to published research, please cite the associated project, manuscript, or repository describing the alignment methodology and spectral retention-time correction framework.

---

## License

TBD
