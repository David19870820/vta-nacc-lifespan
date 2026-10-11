# Data notes

## What is in this repository

| Item | Included? | Notes |
|---|---|---|
| Analysis code | yes | MIT licence |
| Derived subject-level tables | yes | `data_derived/` |
| Result logs, figures | yes | `results/` |
| Raw imaging data (NIfTI) | **no** | fetched by `01_download_bids.py` |
| C-PAC derivatives | **no** | FD values fetched by `03_fetch_fd.py` |

Total repository size is well under 1 MB. This is deliberate: the imaging data
belong to their original distributor and should not be mirrored.

## Source data

**NKI-Rockland** (Enhanced Nathan Kline Institute - Rockland Sample)

- Raw BIDS:
  `s3://fcp-indi/data/Projects/RocklandSample/RawDataBIDSLatest`
- HTTP equivalent:
  `https://fcp-indi.s3.amazonaws.com/data/Projects/RocklandSample/RawDataBIDSLatest/...`
- C-PAC derivatives:
  `s3://fcp-indi/data/Projects/RocklandSample/Outputs/C-PAC`
- Phenotypic file: `aws_links.csv`

The bucket is open and requires no credentials (`--no-sign-request` if using
the AWS CLI). Download size for the analysed subset is roughly 40-50 GB.

## Required citation

Anyone using these data must cite the dataset paper:

> Nooner, K. B., Colcombe, S. J., Tobe, R. H., Mennes, M., Benedict, M. M.,
> Moreno, A. L., Panek, L. J., Brown, S., Zavitz, S. T., Li, Q., Sikka, S.,
> Gutman, D., Bangaru, S., Schlachter, R. T., Kamiel, S. M., Anwar, A. R.,
> Hinz, C. M., Kaplan, M. S., Rachlin, A. B., ... & Milham, M. P. (2012).
> The NKI-Rockland sample: A model for accelerating the pace of discovery
> science in psychiatry. *Frontiers in Neuroscience, 6*, 152.
> https://doi.org/10.3389/fnins.2012.00152

## Derived data in `data_derived/`

Two CSV files are shipped so that the statistics can be reproduced without
downloading any imaging data.

### `primary_onescan_per_subject.csv` — n = 307

| column | meaning |
|---|---|
| `subject` | NKI-Rockland participant ID (e.g. `sub-A00008326`) |
| `age` | age in years at the scan |
| `sex_M` | 1 = male, 0 = otherwise (derived from the `gender` field) |
| `z` | Fisher-z-transformed VTA-NAcc coupling |
| `r` | the same coupling before the z-transform |
| `session` | session label of the retained scan (e.g. `ses-BAS1`) |
| `mean_fd` | mean framewise displacement in mm; NaN where C-PAC had no FD |

These are **derived scalar measures** — one correlation coefficient, one age,
one sex code and one motion summary per participant. No image, no voxel-level
data, no time series, and no directly identifying information is included.

### `full_all_scans.csv` — n = 417

Same columns, one row per successfully processed scan. `session` and `mean_fd`
are blank/NaN where the scan could not be matched to a C-PAC session.

## On redistribution

The derived tables are provided solely to reproduce the published analyses. If
you intend to redistribute them further, or to combine them with other
NKI-Rockland variables, please consult the current terms of use published by
the FCP-INDI and, where required, obtain access through the original channel.

The authors of this repository are not affiliated with the NKI-Rockland study
and make no claim over the underlying data.

## Head-motion source

FD values come from the C-PAC derivatives at

```
Outputs/C-PAC/<sub>/output/pipeline_analysis_nuisance/
    <sub>_<ses>/frame_wise_displacement_power/_scan_rest_acq-1400/FD.1D
```

Note that C-PAC did not process every session; FD was retrievable for 307
subjects and, for the primary scan, for 250 of them.

---

# Appendix: redistribution status of the derived tables

Because this repository ships subject-level derived data, the licensing
position was checked explicitly. Summary of what was found (all verified
against primary sources):

**NKI-Rockland has a two-tier structure.**

- *Open tier* — the neuroimaging release (`RawDataBIDSLatest`, which includes
  age, sex and handedness) plus physiological and motion data. The official
  page states it is "Immediately Available ... No Data Usage Agreement Needed"
  (<https://rocklandsample.org/for-researchers/step-2-data-usage-agreement>).
- *DUA tier* — the full phenotypic release (item-level psychiatric, cognitive
  and behavioural measures) requires a notarised Data Usage Agreement.

**Five of the seven columns shipped here are already published per subject by
NKI-RS itself, openly and anonymously.** The dataset's own motion-statistics
release is a per-subject CSV whose header includes a column literally named
`MeanFD`:

> `Subject,Scan,Mean_Relative_RMS_Displacement,...,MeanFD,NumFD_greater_than_0.20,...`
> <http://fcon_1000.projects.nitrc.org/indi/enhanced/neurodata.html>

Per-subject `age`, `sex`, `handedness` and session labels are likewise public
(`participants.tsv`, `aws_links.csv`, `MIVAC_sampledata.csv` — the last of which
includes race, BMI and education, i.e. considerably more than is shipped here).

The only genuinely new quantities in `data_derived/` are `z` and `r`: one
correlation coefficient per participant, computed from the open-tier images.

**Licence of the open tier.** The AWS Open Data registry lists the FCP-INDI
datasets as `ODC-By v1.0` for imaging data
(<https://registry.opendata.aws/fcp-indi/>). ODC-By is an attribution-only
licence that permits redistribution. Caveats: the bucket itself contains no
`LICENSE` file and `dataset_description.json` carries no licence field, so this
is a registry-level declaration rather than an embedded one; and FCP-INDI's own
access page elsewhere refers to "non-commercial purposes", which is not fully
reconciled with ODC-By.

**No prohibition on redistribution was found for the open tier.** Equally, no
explicit statement about *derived* data was found either — the distinction is
simply not addressed in the public documentation.

**Precedent.** Derived NKI-Rockland connectomes at the subject level have been
published on Zenodo under CC BY (Kiar et al., 2020,
<https://zenodo.org/records/4041549>). No case was found of derived
NKI-Rockland data being taken down.

**Field norms.** The Open Brain Consent data user agreement template addresses
derived-data redistribution in a dedicated clause, but marks it *optional* and
explicitly describes it as "contentious"; its default stance is group-level
only. There is no consensus in the field on subject-level derived measures.

## What this means for the tables here

The tables are provided for reproducibility of the published analyses. Please:

1. Cite Nooner et al. (2012) if you use them.
2. Do not relicense them. **Do not apply CC BY or any other licence to the data
   files** — the code is MIT, but the upstream terms for the data are set by
   NKI-RS/FCP-INDI.
3. Treat them as derived from the NKI-Rockland open neuroimaging tier.

## Open question for the original author

The one variable that could change this assessment is whether a Data Usage
Agreement was signed for the *phenotypic* tier, and what that agreement says
about redistribution. The DUA text is not publicly available; anyone relying on
this repository should satisfy themselves on that point.

If a DUA was signed and prohibits redistribution, the derived tables should be
removed and reproduction limited to the pipeline scripts (`code/01`-`04`),
which regenerate them from the public source. Nothing else in this repository
depends on the tables.
