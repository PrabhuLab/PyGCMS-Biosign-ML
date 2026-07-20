# AlignSpectra

## Overview

Gas chromatography-mass spectrometry (GC-MS) data collected on different instruments often exhibit retention-time shifts, even when analyzing the same sample. These shifts can occur when transitioning between older and newer instruments, after routine maintenance, or when a GC column is trimmed, causing peaks to appear at different scan numbers than before. 

This project aims to provide a framework for aligning py-GC-MS datasets collected under different instrument conditions. The pipeline identifies corresponding peaks between reference datasets and models the relationship between their scan positions using robust systematic correction functions (logarithmic, linear, quadratic, or exponential). Once calibrated, the optimal model can shift historical data so that peaks collected on different instruments align perfectly.

The processing pipeline includes baseline correction, denoising, peak detection, calibration of multiple alignment models, and export of corrected datasets. By aligning data collected before and after maintenance events or across multiple instruments, previously collected samples can be incorporated into the same analysis workflow without requiring manual retention-time adjustments.

---

## Features

* Import py-GC-MS data exported from instrument software or load directly from Pandas DataFrames
* Preserve original metadata headers
* Baseline correction using Asymmetric Least Squares (ALS)
* Signal denoising using urQRd
* Automated peak detection
* Calibration of scan shifts between instruments across multiple mathematical models
* Systematic retention-time alignment (Logarithmic, Linear, Quadratic, and Exponential)
* Visualization of raw and processed chromatograms
* Export of corrected datasets

---

## Installation

### Required Dependencies

```bash
pip install numpy pandas scipy matplotlib
```

### Optional Dependency

Preprocessing requires urQRd, which is located in this folder. Please note that this version of urQRd was modified by Miles Walters from the available copy of the original algorithm found [here](https://doi.org/10.1073/pnas.1306700111).   

Without urQRd, alignment calibration and preprocessing functions will not be available.

---

## Data Format

The software operates on exported py-GC-MS files containing scan-indexed spectra. These can be found in the .3D files in the Py-GC-MS outputs. 

Example:

```text
Abundances from D:\[USER]\GCMS\1\data\[PROJECT]\FILENAME.D   Edit C:\[USER]\MSEXE\\export3d.mac if desired to change output.
Mass Range,     50,    700,
Scan Range,      1,   6295,
Scan number is in the first column
Scan,     50,     51,     52,     53,     54,     55,  ...
...
```

The loader automatically:

* Separates metadata from spectral data
* Uses scan number as the index
* Removes unnamed columns

---

## Basic Usage

### Loading a Dataset

```python
from AlignSpectra import AlignSpectra

processor = AlignSpectra(
    file_path="sample.3D"
)
```

The dataset is automatically loaded during initialization. You can also load an existing Pandas DataFrame directly:

```python
processor.load_df(existing_dataframe)
```

---

### Accessing the DataFrame

```python
df = processor.df

print(df.head())
```

---

## Preprocessing

Preprocessing consists of:

1. Summation of spectral intensities
2. Baseline correction
3. urQRd denoising

```python
processed = processor.preprocess()
```

The processed signal is stored as a numpy array and can be accessed later:

```python
processor.processed_data
```

---

## Plotting

### Processed Chromatogram

```python
import matplotlib.pyplot as plt

processor.preprocess()

processor.plot_spectra()

plt.show()
```

### Raw Chromatogram

```python
processor.plot_spectra(
    processed=False
)

plt.show()
```

---

## Peak Detection

Peak detection identifies local maxima in the processed chromatogram.

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
    step=110,
    threshold=18000
)
```

| Parameter | Description                            |
| --------- | -------------------------------------- |
| step      | Window size used during peak searching |
| threshold | Minimum peak intensity                 |

---

## Calibrating Alignment Parameters

Calibration requires paired datasets collected from different instruments or instrument states. The system will automatically test multiple models (log, linear, quadratic, exponential) to find the best fit.

Example:

```python
pairs = [
    (
        "old_machine_run_1.3D",
        "new_machine_run_1.3D"
    ),
    (
        "old_machine_run_2.3D",
        "new_machine_run_2.3D"
    )
]

best_model, best_params, window, accuracy = (
    AlignSpectra.calibrate_shift_parameters(
        pairs,
        verbose=True
    )
)

print("Optimal model:", best_model)
print("Optimal window:", window, "scans")
print("Optimal parameters:", best_params)
print(f"Accuracy: {accuracy*100:.2f}%")
```

Output:

```text
Optimal model: log
Optimal window: 85
Optimal parameters: {'a': 4.83, 'b': 2.14}
Accuracy: 96.80%
```

### Calibration Method

For each paired dataset:

1. Preprocess both chromatograms
2. Detect chromatographic peaks
3. Match corresponding peaks
4. Measure scan offsets
5. Remove outliers
6. Fit multiple predictive models (logarithmic, linear, quadratic, exponential)

The method iterates over peak detection window sizes and model types to determine the optimal mathematical approach for creating the most accurate shift. 

---

## Applying Alignment

After calibration, you can apply the learned model to historical datasets:

```python
processor = AlignSpectra(
    file_path="old_sample.3D",
    is_old=True
)

processor.set_shift_parameters(
    best_model,
    best_params
)
```

Depending on the calibrated model, the alignment function uses one of the following formulas:

* **Logarithmic:** `shifted_scan = scan + (a * np.log(scan) + b)`
* **Linear:** `shifted_scan = scan + (a * scan + b)`
* **Quadratic:** `shifted_scan = scan + (a * scan**2 + b * scan + c)`
* **Exponential:** `shifted_scan = scan + (a * np.exp(b * scan) + c)`

---

## Exporting Aligned Data

```python
processor.save_shifted(
    "aligned/aligned_sample.3D"
)
```

The exported file:

* Preserves metadata
* Applies scan correction based on the chosen model
* Writes a new .3D file

---

## Complete Example

```python
from AlignSpectra import AlignSpectra

training_pairs = [
    (
        "old_reference.3D",
        "new_reference.3D"
    )
]

best_model, best_params, _, accuracy = (
    AlignSpectra.calibrate_shift_parameters(
        training_pairs,
        verbose=True
    )
)

print(
    f"Calibration accuracy: {accuracy:.2%}"
)

sample = AlignSpectra(
    file_path="old_unknown.3D",
    is_old=True
)

sample.set_shift_parameters(
    best_model,
    best_params
)

sample.preprocess()

sample.plot_spectra()

sample.save_shifted(
    "aligned_unknown.3D"
)
```

---

## Constructor Parameters

```python
AlignSpectra(
    file_path=None,
    is_old=True,
    a=None,
    b=None,
    columns=None,
    df=None,
    shift_model=None,
    shift_params=None
)
```

| Parameter    | Description                                                     |
| ------------ | --------------------------------------------------------------- |
| file_path    | Path to spectral file                                           |
| is_old       | Whether alignment should be applied                             |
| a            | *(Deprecated)* Logarithmic coefficient                          |
| b            | *(Deprecated)* Constant offset                                  |
| columns      | Spectral column range used in analysis (50-700 if blank)        |
| df           | Optional existing Pandas DataFrame to load directly             |
| shift_model  | String name of shift model (`'log'`, `'linear'`, `'quadratic'`, `'exponential'`) |
| shift_params | Dictionary of coefficients for the model (e.g., `{'a':1, 'b':2}`) |

---

## Processing Workflow

```text
Raw py-GC-MS Data --> Intensity Summation --> Baseline Correction --> urQRd Denoising --> Peak Detection --> Calibration --> Model Selection & Fitting --> Peak Alignment --> Corrected Dataset
```

---

## Error Handling

### Missing File

```python
ValueError:
No file path provided
```

### Invalid File Format

```python
ValueError:
File format error: 'Scan,' header not found
```

### Missing urQRd

```python
ImportError:
urQRd package is required for preprocessing
```

### Missing Shift Parameters

```python
ValueError:
Shift parameters have not been set.
```

### Unknown Model Selected

```python
ValueError:
Model must be one of ['log', 'linear', 'quadratic', 'exponential']
```
*(or)*
```python
ValueError:
Unknown shift model: '...'
```

---

## Citation

If this software contributes to published research, please cite the associated project, manuscript, or repository describing the alignment methodology and py-GC-MS retention-time correction framework.

---

## License

TBD
