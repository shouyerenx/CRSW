# CRSW-PL

This repository contains the MATLAB implementation of CRSW-PL for
interferometric phase linking of long SAR time series.

The method is described in the manuscript *Interferometric phase linking for
long time series based on coherence magnitude matrix regularization and
overlap fusion*.

## Requirements

- MATLAB R2024b (used for validation) and Image Processing Toolbox.
- Parallel Computing Toolbox for parallel processing (`cfg.workers > 0`).
- Statistics and Machine Learning Toolbox when `cfg.shp_alpha` is neither
  0.05 nor 0.01.

## DGB dataset

The DGB dataset is available on [Zenodo (review link)](https://zenodo.org/records/24298081?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjFmNWM2ZjkxLTk4ODUtNDNjMi05OTJlLTcxOGI0ZjU5Y2ZiYiIsImRhdGEiOnt9LCJyYW5kb20iOiIyOTk0ZGUxZDRjNTEyYjEwM2YyNTJhNjk2ZWQ5ZDQ0MyJ9.XRIcd8HaRuuO9R9KTWFioVTjZPyKanfAvutZjjezSvyOd5q8yTZhpZrpV_fspAmVzHClBJitHAopxqOK6kxv9A).
Download all seven ZIP files and extract each one into the same `DGB` folder
inside the project directory. The resulting layout should be:

```text
CRSW/
  DGB/
    slc/
    rmli/
    diff/
```

For this dataset, set the following values in `CRSW_settings.m`:

```matlab
cfg.input_dir = 'DGB';
cfg.master_date = '20170224';
cfg.nlines = 681;
cfg.slc_format = 'cpxint16';
```

## Input

Set `cfg.input_dir` to a directory containing:

```text
input_dir/
  YYYYMMDD.rslc           complex int16 or complex float32
  YYYYMMDD.rmli           float32 intensity
  MASTER_YYYYMMDD.diff    complex float32
```
The files can also be placed in `slc`, `rmli`, and `diff` subdirectories.

## Run

Set `input_dir`, `output_dir`, `master_date` (in `YYYYMMDD` format), `nlines`,
and `slc_format` in `CRSW_settings.m`. Set `workers` to 0 for serial processing
or to the desired number of workers. In MATLAB, change to the CRSW directory
and run:

```matlab
CRSW
```

The optional crop is `[first row, last row, first column, last column]` with
one-based inclusive indices.

## Output

`output_dir/CRSW_result.mat` contains:

- `patch_ph`: wrapped phase in radians, indexed by row, column, and date, with
  the first acquisition as the phase reference.
- `dates`: acquisition dates.
- `master_date`: input master date.
- `shp_counts`: selected-pixel counts for each temporal window.
- `cfg`: run settings.

## Display an interferogram

Run `showIntf` from the CRSW directory to display the Interferogram:

```matlab
showIntf
```

The script reads `./CRSW_output/CRSW_result.mat`. Edit the path and the
one-based date indices `m` and `n`  in `showIntf.m` as
needed.
