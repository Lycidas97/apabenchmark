# Shape Core Params Plot Debug Guide

This guide focuses on the shape-core boxplot script:
- `scripts/06_plotting/peak_model_params/scripts/25_panel_shape_core_params.py`

## 1) Quick Run (default output)

```bash
python3 scripts/06_plotting/peak_model_params/scripts/25_panel_shape_core_params.py
```

Default output:
- `scripts/06_plotting/peak_model_params/figures/shape_core_params.pdf`

## 2) Trial Run with custom figure size

Use this when tuning layout repeatedly:

```bash
python3 scripts/06_plotting/peak_model_params/scripts/25_panel_shape_core_params.py \
  --width-mm 170 \
  --height-mm 45 \
  --wspace 0.12 \
  --output scripts/06_plotting/peak_model_params/figures/shape_core_params_trial_v1.pdf
```

Key knobs:
- `--fig-preset`: use style presets (default: `mid_wide`)
- `--width-mm` + `--height-mm`: manual size override (must be provided together)
- `--wspace`: horizontal spacing between subplots
- `--output`: write to a dedicated trial file for comparison

## 3) Full topic rerun

```bash
python3 scripts/06_plotting/peak_model_params/scripts/run_all.py
```

This regenerates all panels in the peak-model-params topic.

## 4) Parameters currently plotted

The script plots exactly these 3 metrics:
- `shape_per_peak_upstream_asymmetry_log2_mean`
- `shape_per_peak_effective_width_q10_q90_mean`
- `shape_per_peak_average_distance_mean`

To change metric selection, edit `PARAMETERS` in:
- `scripts/06_plotting/peak_model_params/scripts/25_panel_shape_core_params.py`

## 5) Style/alignment notes

Current micro-style behavior (aligned with reference notebooks and topic scripts):
- inward ticks on both axes
- minor ticks on numeric x-axis
- extra plot margins (`ax.margins(y=0.06, x=0.07)`) to avoid overlap between inward ticks and box elements

If box/tick overlap appears again, tune margins first before changing box width.
