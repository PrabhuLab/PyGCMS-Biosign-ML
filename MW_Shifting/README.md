# AlignSpectra

## Overview

Gas chromatography-mass spectrometry (GC-MS) data collected on different instruments often exhibit retention-time shifts, even when analyzing the same sample. These shifts can occur when transitioning between older and newer instruments, after routine maintenance, or when a GC column is trimmed, causing peaks to appear at different scan numbers than they did previously. As a result, direct comparison between datasets becomes difficult because corresponding compounds no longer align.

This project aims to provide a framework for aligning py-GC-MS datasets collected under different instrument conditions. The software identifies corresponding peaks between reference datasets and models the relationship between their scan positions using a logarithmic correction function. Once calibrated, the model can shift historical data so that peaks collected on different instruments align with one another.

The processing pipeline includes baseline correction, denoising, peak detection, calibration of alignment parameters, and export of corrected datasets. By aligning data collected before and after maintenance events or across multiple instruments, previously collected samples can be incorporated into the same analysis workflow without requiring manual adjustment of retention times.

Understanding how retention times shift between instruments enables the integration of data collected over long periods of time and across different laboratories, reducing instrument-specific variability and improving the consistency of downstream analyses.

---

## Features

* Import py-GC-MS data exported from instrument software
* Preserve original metadata headers
* Baseline correction using Asymmetric Least Squares (ALS)
* Signal denoising using urQRd
* Automated peak detection
* Calibration of scan shifts between instruments
* Logarithmic retention-time alignment
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
Abundances from D:\[USERNAME]\GCMS\1\data\[PROJECT]\FILENAME.D   Edit C:\USERNAME\MSEXE\\export3d.mac if desired to change output.
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
* Restricts analysis to the selected mass range

Default spectral range:

```python
50-700
```

---

## Basic Usage

### Loading a Dataset

```python
from align_spectra import AlignSpectra

processor = AlignSpectra(
    file_path="sample.3D"
)
```

The dataset is automatically loaded during initialization.

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

The processed signal is stored as a DataFrame:

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
    step=80,
    threshold=15000
)
```

| Parameter | Description                            |
| --------- | -------------------------------------- |
| step      | Window size used during peak searching |
| threshold | Minimum peak intensity                 |

---

## Calibrating Alignment Parameters

Calibration requires paired datasets collected from different instruments or instrument states.

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

a, b, window, accuracy = (
    AlignSpectra.calibrate_shift_parameters(
        pairs,
        verbose=True
    )
)

print("Optimal window:", window, "scans")
print(f"Optimal model: {a.round(2)}*ln(x) + {b.round(2)}"
print(f"Accuracy: {accuracy}%")
```

Output:

```text
Optimal window: 85 scans
Optimal model: 4.83*ln(x) + 2.14
Accuracy: 96.8%
```

### Calibration Method

For each paired dataset:

1. Preprocess both chromatograms
2. Detect chromatographic peaks
3. Match corresponding peaks
4. Measure scan offsets
5. Remove outliers
6. Fit a logarithmic model

The resulting model has the form:

offset = a * ln(x) + b

where:

* x is the scan number
* a and b are learned parameters

---

## Applying Alignment

After calibration:

```python
processor = AlignSpectra(
    file_path="old_sample.3D",
    is_old=True
)

processor.set_shift_parameters(
    a,
    b
)
```

The alignment function is:

```python
shifted_scan = scan + (
    a * np.log(scan) + b
)
```

---

## Exporting Aligned Data

```python
processor.save_shifted(
    "aligned/aligned_sample.3D"
)
```

The exported file:

* Preserves metadata
* Applies scan correction
* Writes a new CSV file

---

## Complete Example

```python
from align_spectra import AlignSpectra

training_pairs = [
    (
        "old_reference.3D",
        "new_reference.3D"
    )
]

a, b, _, accuracy = (
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
    a,
    b
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
    columns=None
)
```

| Parameter | Description                                              |
| --------- | -------------------------------------------------------- |
| file_path | Path to spectral file                                    |
| is_old    | Whether alignment should be applied                      |
| a         | Logarithmic coefficient                                  |
| b         | Constant offset                                          |
| columns   | Spectral column range used in analysis (50-700 if blank) |

---

## Processing Workflow

```text
Raw py-GC-MS Data --> Intensity Summation --> Baseline Correction --> urQRd Denoising --> Peak Detection --> Calibration --> Logarithmic Shift Model --> Peak Alignment --> Corrected Dataset
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
File format error:
'Scan,' header not found
```

### Missing urQRd

```python
ImportError:
urQRd package is required for preprocessing
```

### Missing Shift Parameters

```python
ValueError:
No valid "a" value
```

```python
ValueError:
No valid "b" value
```

---

## Citation

If this software contributes to published research, please cite the associated project, manuscript, or repository describing the alignment methodology and py-GC-MS retention-time correction framework.

---

## License

TBD
