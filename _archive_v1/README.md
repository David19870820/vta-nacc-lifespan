# Archived: first-version (v1) artefacts

These files were produced by the **first version** of the analysis, before
head motion was obtained and before the equivalence / reliability analyses
were added. They are kept for provenance and are **not** part of the
current manuscript.

Do not cite these values. In particular `vta_nacc_stats.json` carries the
note `"note_head_motion": "NOT controlled (no confound in raw BIDS)"`,
which was true of v1 but is **false** of the current analysis: FD values
were subsequently obtained from the C-PAC derivatives and are controlled
throughout Sections 3.4 and 3.7. See `../MANUSCRIPT.md` and
`../DATA_NOTES.md`.

Current equivalents:

| archived file | current equivalent |
|---|---|
| `vta_nacc_noage.csv` | `data_derived/full_all_scans.csv` |
| `vta_nacc_stats.json` / `.txt` | `results/quick_reproduce.txt` |
| `vta_nacc_age_scatter.png` | `results/Fig1_age_vs_coupling.png` |
