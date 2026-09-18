# Input data format

The ready-to-run file [`Example data.xlsx`](Example%20data.xlsx) demonstrates the exact input layout accepted by `run_analysis.m`.

## Layout

The first worksheet must contain a numeric matrix with no header row and no blank cells.

| Position | Meaning | Required values |
| --- | --- | --- |
| Column 1 | Trial coordinate | Finite, unique, strictly increasing numbers |
| Columns 2 onward | Freezing for one animal per column | Finite percentages from 0 to 100 |
| Rows | Measured trials | At least five rows for the manuscript model comparison |

The supplied example has five rows and 23 columns: one trial-coordinate column followed by 22 animal columns. The coordinates are `0, 1, 2, 3, 4`, corresponding to displayed trial numbers 1-5. The code linearly interpolates these measurements to a spacing of 0.2 before fitting. The fitted tdRL value `mu_trial` is therefore reported as `mu + 1`.

Minimal schematic example:

```text
0   animal_01_trial_1   animal_02_trial_1   ...
1   animal_01_trial_2   animal_02_trial_2   ...
2   animal_01_trial_3   animal_02_trial_3   ...
3   animal_01_trial_4   animal_02_trial_4   ...
4   animal_01_trial_5   animal_02_trial_5   ...
```

## Validation

`prepare_behavior_data.m` stops with a specific error if the input file is missing, contains headers/text/missing values, has invalid dimensions, contains out-of-range freezing percentages, or uses unordered trial coordinates. This prevents a malformed file from silently entering the fit.

To analyze a different group, prepare a file with the same layout and run:

```matlab
run_analysis('C:\path\to\group_data.xlsx', 'C:\path\to\group_results')
```

Fit fast- and slow-learning groups separately, using a separate output directory for each group.
