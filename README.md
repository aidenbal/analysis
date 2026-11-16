# White Matter Tractometry, Stress, and Cardiovascular Health

Diffusion MRI tractometry (PyAFQ) and Elastic Net prediction applied to the
LEMON dataset, testing whether white matter microstructure is associated with
perceived stress, chronic stress, and mean arterial pressure (MAP).

## Research question

Do white matter pathways implicated in stress regulation show localized associations
with stress-related measures, and can they predict these outcomes at the individual
level beyond age and sex?

## Data

This project uses the **LEMON** dataset (Leipzig Study for Mind-Body-Emotion
Interactions). Participant-level data are **not** included in this repository;
they can be obtained from the original data providers under their terms of use.

| Outcome | Type |
|---|---|
| PSQ (Perceived Stress Questionnaire, overall score) | Self-reported perceived stress |
| TICS (Trier Inventory for Chronic Stress, screening scale) | Self-reported chronic stress |
| MAP (mean arterial pressure) | Physiological, stress-related outcome (not a direct stress measure) |

MAP = diastolic + (systolic − diastolic) / 3, using the mean of all available readings.

## Analyses

1. **Preprocessing (QSIPrep):** denoising, susceptibility distortion correction,
   eddy-current and motion correction, b-vector rotation, registration to anatomical space.
2. **Tractometry (PyAFQ):** bundle segmentation and FA profiles at 100 nodes per tract.
3. **Analysis 1 – Whole-brain node-wise screening:** at every node of every tract with
   at least 10 participants, `FA ~ outcome + Age + Sex`, with Benjamini-Hochberg FDR
   correction across nodes.
4. **Analysis 2 – ROI node-wise analysis:** three tract systems implicated in stress
   regulation: corpus callosum (CCMid), bilateral uncinate fasciculus (UF_L, UF_R), and
   bilateral superior longitudinal fasciculus I–III (six tracts). Clusters are runs of at
   least 3 consecutive nodes with p < .05.
5. **Analysis 3 – Elastic Net prediction (glmnet):** mean FA of the nine ROI tracts plus
   age and sex (alpha = 0.5). Performance is estimated with **nested** 10-fold
   cross-validation (lambda tuned inside each training split), repeated 20 times.
6. **Analysis 4 – Baseline comparison:** the same cross-validation for an age + sex
   model, on the same participants and folds, to test whether white matter adds
   predictive value.
7. **Analysis 5 – Ordinary regression and joint test (MAP):** coefficients with p-values
   and an F-test of all tracts together beyond age and sex.

## Key results

**Node-wise screening and ROI analysis.** Exploratory associations appeared in CCMid,
bilateral UF and the SLF (for example, CCMid with PSQ and TICS, UF with all three
outcomes, left and right SLF II with MAP, right SLF II with PSQ). Only MAP nodes in Callosum Anterior Frontal (N = 11; not interpreted) survived FDR correction. Full cluster
lists are in `results/roi_clusters.csv`.

**Prediction (nested cross-validated, mean of 20 repeats).**

| Outcome | Model | N | CV R² | CV RMSE |
|---|---|---|---|---|
| MAP | Age + sex only | 49 | 0.43 | 10.7 mmHg |
| MAP | 9 tracts + age + sex | 49 | 0.35 | 11.4 mmHg |
| PSQ | Age + sex only | 45 | −0.04 | 17.5 |
| PSQ | 9 tracts + age + sex | 45 | 0.14 | 15.9 |
| TICS | Age + sex only | 45 | −0.07 | 8.9 |
| TICS | 9 tracts + age + sex | 45 | −0.02 | 8.7 |

MAP was predictable, but through age; adding white matter did not improve prediction
beyond age and sex. For PSQ, the white matter model explained about 14% of variance
where age and sex explained none, an exploratory result given the small sample.
TICS could not be predicted.

**Take-home.** White matter in the corpus callosum, uncinate fasciculus and SLF showed
exploratory node-wise associations with stress-related measures. These tracts did not
improve MAP prediction beyond age and sex, but showed modest predictive value for
perceived stress (PSQ) that warrants replication in a larger sample.

## Sample sizes

**Node-wise analyses (Analyses 1–2).** Participants per tract. Sample sizes vary because
PyAFQ did not reconstruct every tract in every participant. Tracts found in fewer than 10
participants (AF_L, AF_R, CST_L, CST_R, CC_ForcepsMajor, Callosum Motor, Callosum Superior
Frontal) were not analysed.

| Tract | MAP | PSQ | TICS | Note |
|---|---|---|---|---|
| C_L | 212 | 205 | 205 | |
| Right Superior Longitudinal III | 211 | 205 | 205 | ROI tract |
| Right Superior Longitudinal II | 194 | 188 | 188 | ROI tract |
| Left Superior Longitudinal III | 191 | 186 | 186 | ROI tract |
| Left Superior Longitudinal II | 179 | 173 | 173 | ROI tract |
| MCP | 169 | 163 | 163 | |
| UF_R | 164 | 158 | 158 | ROI tract |
| C_R | 162 | 156 | 156 | |
| UF_L | 151 | 146 | 146 | ROI tract |
| CCMid | 148 | 143 | 143 | ROI tract |
| Right Superior Longitudinal I | 145 | 138 | 138 | ROI tract |
| Left Superior Longitudinal I | 142 | 137 | 137 | ROI tract |
| CC_ForcepsMinor | 102 | 97 | 97 | |
| F_R | 23 | 23 | 23 | Too few participants; not interpreted |
| F_L | 22 | 22 | 22 | Too few participants; not interpreted |
| Callosum Anterior Frontal | 11 | 11 | 11 | Too few participants; not interpreted |

**Prediction analyses (Analyses 3–5).** These required complete data for all nine ROI
tracts: N = 49 for MAP and N = 45 for PSQ and TICS.

## Notes and limitations

- An earlier version of this analysis used four tracts (CCMid, UF_L, UF_R, right SLF II)
  and reported R² ≈ 0.515 for MAP. That value was an in-sample fit; nested
  cross-validation gave 0.41, and age + sex alone gave 0.43.
- Requiring all nine tracts reduces the prediction sample to 45–49 participants, so
  estimates are imprecise.
- Results are exploratory and need replication in an independent sample.

## Repository structure

```
analysis_code.R   Full analysis: data merge, node-wise regressions, ROI clusters,
                  nested-CV Elastic Net, baseline comparison, regression, figures
results/          Summary tables (CSV)
figures/          Figures (PNG)
```

## How to reproduce

1. Obtain the LEMON data and run QSIPrep and PyAFQ to produce tract profiles.
2. Place the behavioural/BP CSVs and the two PyAFQ profile CSVs in the project folder
   (file names are listed at the top of `analysis_code.R`).
3. Open `analysis.Rproj` in RStudio and run `source("analysis_code.R")`.
   Required packages (dplyr, tidyr, ggplot2, glmnet) install automatically.

Software versions are recorded in `results/session_info.txt`.

## Author

Aiden Balandan
