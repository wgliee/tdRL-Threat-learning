# R-W and trial-wise discount reinforcement-learning models

This repository fits and compares two models of freezing across threat-learning trials:

1. the standard Rescorla-Wagner (R-W) model; and
2. the trial-wise discount reinforcement-learning (tdRL) model described in the accompanying manuscript.

The repository includes a ready-to-run input file, automatic input validation, an end-to-end example, and machine-readable output files. No manual conversion of the supplied example data is required.

## Requirements

The analysis was developed and tested for **MATLAB R2023b** and uses:

- MATLAB (`ode45`, data import, tables, and plotting);
- Optimization Toolbox (`optimvar`, `optimproblem`, `fcn2optimexpr`, and `solve`); and
- Statistics and Machine Learning Toolbox (`lillietest`, `ttest`, and `signrank`).

Run `check_dependencies` in MATLAB to verify that the required functions are available.

## Quick start

1. Download or clone the complete repository. Keep `Example data.xlsx`, the `.m` files, and the `function` folder together.
2. Open MATLAB and set the current folder to this repository root.
3. Run:

```matlab
run_example
```

The command checks dependencies, validates the example input, fits R-W, fits tdRL, compares their animal-level errors, and writes all outputs to `results/`.

For a custom data file and output directory, use:

```matlab
run_analysis('C:\path\to\data.xlsx', 'C:\path\to\results')
```

Paths may point to `.xlsx`, `.xls`, `.csv`, or numeric text files supported by MATLAB `readmatrix`. See [DATA_FORMAT.md](DATA_FORMAT.md) for the exact row/column structure and validation rules.

## Repository contents

| File or folder | Purpose |
| --- | --- |
| `Example data.xlsx` | Ready-to-run example: five measured trials from 22 animals |
| `run_example.m` | Runs the complete example with no arguments |
| `run_analysis.m` | Runs the complete workflow for a specified input and output path |
| `check_dependencies.m` | Checks MATLAB release and required functions |
| `RW_G.m` | Fits the R-W model and exports its predictions and parameters |
| `tdRL_G.m` | Fits the tdRL model and exports its predictions and parameters |
| `compare_RW_vs_tdRL_G.m` | Calculates animal-level errors, selects the paired test, and exports results |
| `plot_pdf_G.m` | Plots the normalized Gaussian discount function used by tdRL |
| `function/` | ODE right-hand sides, multi-animal simulators, and data helpers |
| `DATA_FORMAT.md` | Input data dictionary, example layout, and custom-data instructions |

## Models and variables

For R-W, the value estimate evolves as

```text
dV/dt = alpha * (lambda - V)
```

For tdRL, learning is reduced by a normalized Gaussian trial-wise discount function:

```text
phi(t) = exp(-(t - mu)^2 / (2 * sigma^2))
dV/dt = (1 - phi(t)) * alpha * (lambda - V)
```

| Variable | Description | Fitting bounds |
| --- | --- | --- |
| `V` | Predicted freezing (%) | Model state |
| `t` | Trial coordinate; example data use 0-4 for displayed trials 1-5 | Input/model grid |
| `lambda` | Asymptotic freezing level (%) | 0-100 |
| `alpha` | Learning rate | 0-1 |
| `mu` | Center of the tdRL discount function on model coordinates | 0-6 |
| `mu_trial` | Center reported on displayed trial numbers (`mu + 1`) | Derived value |
| `sigma` | Positive width of the tdRL discount function | greater than 0 and at most 20 |

The initial parameter vectors are `[50, 0.1]` for R-W and `[50, 0.1, 0.1, 0.4]` for tdRL, matching the manuscript Methods.

## Analysis workflow

The five measured trial values are linearly interpolated at intervals of 0.2. Each animal's first observed value initializes its ODE trajectory. Parameters are shared within the analyzed group and are estimated by minimizing the sum of squared differences between the predicted and observed group means. Individual predicted trajectories are then generated with those fitted parameters.

Model performance is compared at the first five measured trials. Freezing percentages are converted to proportions, and SSE, MSE, and RMSE are calculated separately for every animal. The script applies a Lilliefors test to each model's animal-level RMSE distribution. If both distributions pass, it uses a paired t-test; otherwise, it uses a Wilcoxon signed-rank test. This implements the decision rule described in the manuscript.

The manuscript analyzes fast- and slow-learning groups separately. To reproduce those group-specific results, call `run_analysis` once for each group and use a different output directory for each call.

## Outputs

The generated `results/` directory contains:

- `RW_fit_parameters.csv` and `tdRL_fit_parameters.csv`: fitted parameters, objective SSE, and solver exit flag;
- `RW_predictions.csv` and `tdRL_predictions.csv`: interpolated trial coordinates and one predicted column per animal;
- `observed_interpolated.csv`: linearly interpolated observed trajectories;
- `RW_fit_results.mat` and `tdRL_fit_results.mat`: complete fit inputs and outputs;
- `model_comparison_per_animal.csv`: each animal's SSE, MSE, and RMSE for both models;
- `model_comparison_summary.csv`: mean and SD of RMSE with the animal count;
- `model_comparison_normality.csv`: Lilliefors decisions and p-values;
- `model_comparison_test.csv`: selected paired test, statistic, degrees of freedom where applicable, p-value, confidence interval where applicable, and alpha;
- `model_comparison_results.mat`: complete comparison variables; and
- PNG figures for both fits, residuals, residual distributions, and cumulative error.

The workflow does not use random sampling, so repeated runs with the same MATLAB/toolbox versions and identical input should be deterministic apart from possible platform-level numerical tolerances.

## Troubleshooting

- **Input file not found:** confirm that the full path passed to `run_analysis` exists. `run_example` expects `Example data.xlsx` in the repository root.
- **Non-numeric or missing data:** remove headers, labels, blank cells, `NaN`, and `Inf`. Every cell in the input range must be numeric.
- **Missing MATLAB function:** run `check_dependencies` and install or enable the toolbox named in the error.
- **Comparison result files missing:** run `run_example`, or run `RW_G.m` and `tdRL_G.m` before `compare_RW_vs_tdRL_G.m`.
- **Different experimental groups:** use separate input files and output directories. Do not combine fast- and slow-learning animals in one fit if reproducing the manuscript analysis.
