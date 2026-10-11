# Resting-State Functional Coupling Between the Ventral Tegmental Area and Nucleus Accumbens Shows No Detectable Age Effect Across the Lifespan: A Replication in the NKI-Rockland Sample

Xiaowei Hu, Independent Researcher

## Abstract

**Background.** The mesolimbic dopamine system is central to motivated behavior. Murty et al. (2018) reported that functional coupling between the ventral tegmental area (VTA) and the nucleus accumbens (NAcc) varies with age in a motivational context but **not** in a resting-state context. Their sample spanned 180 participants aged 10–30 years in a longitudinal design. A null effect at rest has since been noted in two further developmental samples (Tomasi & Volkow, 2014; Morales et al., 2024), but it has never been formally replicated, never been tested beyond roughly 30 years of age, and never been examined against head motion, a confound that increases with age and can reduce apparent connectivity.

**Methods.** We analysed resting-state fMRI (TR = 1400 ms) from 307 unique participants (417 scans; primary analysis used one scan per participant), age 9–83 years, from the publicly available NKI-Rockland sample. VTA and bilateral NAcc were defined as 4-mm spheres in MNI space (VTA: 0, −16, −14; NAcc: ±10, 12, −8). Coupling was the Fisher-z-transformed Pearson correlation between VTA and mean NAcc time series after detrending, standardizing, and band-pass filtering (0.01–0.08 Hz). Age effects were tested with linear, inverse-age (1/age, replicating Murty et al.), and multiple-regression models controlling for sex. Head motion (framewise displacement, FD) was obtained from C-PAC derivatives for the same acquisition and controlled in three ways: as a covariate, by excluding high-motion participants, and by diagnostic tests of the motion–age and motion–coupling relationships.

**Results.** There was no age effect on VTA–NAcc resting-state coupling. Linear correlation: r = −0.042, p = 0.46 (n = 307). Inverse-age model with sex: β = −0.703, p = 0.34. Trimming the youngest and oldest 5% and 10% left p > 0.40. Across 15 model specifications the partial correlation of age with coupling never exceeded |r| = 0.10 (Figure 3). Head motion increased modestly with age (r = +0.18, p = 0.005) but was unrelated to coupling (r = +0.06, p = 0.38, n = 250). Adding mean FD as a covariate left the age effect non-significant (p = 0.82; inverse-age p = 0.64), as did excluding high-motion participants at five thresholds (all p > 0.32), and the full-scan analysis (304 scans, 269 participants; p = 0.75). A linear mixed model with a random intercept per participant gave the same result (age p = 0.75). Equivalence testing (TOST) excludes age effects of |r| ≥ 0.14; smaller effects remain compatible with the data. A Bayes factor of BF01 = 8.5–16.7 constitutes moderate-to-strong evidence for the absence of an age effect. A secondary finding was higher coupling in male than female participants (β = 0.073, p = 0.019).

**Conclusion.** VTA–NAcc resting-state functional coupling shows no detectable age effect from childhood to old age in a large lifespan sample, formally replicating and extending Murty et al. (2018). The null is not attributable to head motion. The age-dependence of VTA–NAcc coupling reported in motivational contexts therefore does not reflect a tonic property of this specific mesolimbic edge. We emphasise that this age-invariance is edge-specific: within the same dataset, dorsal striatal connectivity topography declines with age, and other VTA resting-state connections are known to change.

**Robustness.** Recomputing coupling from fully preprocessed C-PAC derivatives (nuisance-regressed, motion-corrected) in 178 participants left the age effect non-significant (r = +0.045, p = 0.55); the two pipelines' coupling estimates correlated at r = 0.28, and test–retest reliability was r = 0.28 (ICC = 0.26) across 40 same-type session pairs and r = 0.11 (ICC = 0.12) across all 182 session pairs. Across twelve ROI definitions (radii 2–6 mm, ±2 mm seed displacements, an alternative NAcc target) no zero-order age correlation was significant (all p > 0.26).

**Limitations.** The main analysis used raw BIDS images without nuisance regression or volume scrubbing, with FD modelled as a covariate rather than removed from the estimates. FD was unavailable for 57 of 307 participants, and preprocessed derivatives for 129 of 307. Single-scan reliability is fair at best (ICC = 0.26 for same-type pairs, 0.12 across all pairs). After correcting for attenuation, the equivalence claim excludes only large true effects (|ρ| ≥ 0.27), not medium ones.

**Keywords:** ventral tegmental area; nucleus accumbens; resting-state functional connectivity; lifespan development; replication; head motion

## 1. Introduction

The mesolimbic dopamine system, comprising the ventral tegmental area (VTA) and its projections to the nucleus accumbens (NAcc), is a core substrate of motivation and reward processing. A key developmental question is whether this system's engagement changes across the lifespan. Adolescence is characterized by heightened motivation and increased sensitivity to reward, which has been attributed to developmental changes in mesolimbic dopamine signaling. Because dopaminergic markers such as dopamine transporter and D2 receptor availability decline from early adulthood onward, the question extends naturally into middle and old age.

### 1.1 The resting-state null and what remains untested

Murty et al. (2018) addressed this question directly using functional connectivity. In a longitudinal sample of 180 participants (age 10.1–30.3 years), they examined VTA–NAcc functional coupling in two contexts: a resting-state scan and a motivated (rewarded antisaccade) task. Their central finding was a **context-dependent** developmental trajectory: VTA–NAcc coupling was enhanced in pre-adolescence/early adolescence and decreased into adulthood in the *motivational* context (p = 0.0004), but showed **no significant age effect in the resting-state context** (p = 0.56). This dissociation suggests that age-related differences in mesolimbic engagement are state-dependent. These differences emerge when the system is recruited by goal-relevant cues, but are absent in its tonic, task-free state.

The resting-state null has not gone unremarked, and we wish to be explicit that the present study does not report it for the first time. Tomasi and Volkow (2014), comparing 402 children/adolescents with 704 adults, reported age-related *increases* in VTA connectivity with several limbic regions, but their VTA–NAcc contrast did not reach significance (right NAcc: not significant for the VTA seed, though significant for the substantia nigra seed). Morales et al. (2024), in a nine-wave longitudinal cohort (N = 796, 4,025 scans, baseline age 12–21), reported the VTA–NAcc age term as β = −0.006, p = 0.83, and described it in their introduction as showing “no significant associations between VTA-NAcc RSFC from adolescence to adulthood.”

What has not been done is the following. First, **no study has formally replicated the null**, in the sense of pre-specifying it as the target of inference rather than reporting it as a background or incidental finding; a search of the citing literature of Murty et al. (2018) returns no such attempt. Second, **no study has tested it beyond approximately 30 years of age.** Both the original sample and the two subsequent developmental cohorts stop at or before the transition to middle adulthood, precisely where dopaminergic decline becomes pronounced. Third, and most importantly for a null result, **motion has not been addressed.** Head motion increases with age and systematically reduces apparent functional connectivity; an uncorrected null is therefore open to the objection that it reflects motion-driven loss of sensitivity rather than a genuinely age-invariant coupling.

The present study targets exactly these three gaps.

### 1.2 Why stability would be informative, and why it is not a claim about “the mesolimbic system” in general

If resting-state VTA–NAcc coupling is stable well into old age, while dopaminergic molecular markers decline over the same span, the two facts are in tension with the default expectation in the cognitive ageing literature, where the prevailing framework is **dedifferentiation**: age-related loss of functional specificity, and within the striatum in particular. Rieckmann et al. (2018) showed that striatal dopamine transporter density and caudate functional connectivity dedifferentiation jointly predict memory decline in normal ageing, that is, molecular decline and connectivity change are typically observed together.

Two findings make it essential to state the claim narrowly rather than as a general statement about mesolimbic or striatal intrinsic connectivity:

**Dorsal striatum.** O'Rawe and Leung (2022), using the **same NKI-Rockland data and the same TR = 1400 ms acquisition** analysed here, found robust age-related degradation of caudate connectivity topography (medial–lateral gradient: r = −0.34, p < 0.001 for the 189-participant cohort; r = −0.37, p < 0.001 for the 299-participant enhanced cohort).

**Other VTA edges.** Peterson et al. (2017), in 250 adults aged 18–49, reported significant age effects on VTA resting-state connectivity with superior temporal and parahippocampal regions and the cerebellum (increasing) and with central sulcus and postcentral gyrus (decreasing).

We therefore frame the hypothesis at the level of a single edge: **the tonic VTA–NAcc coupling is developmentally conserved, whereas age effects on the mesolimbic system are expressed in task-evoked recruitment and in other, predominantly dorsal, connections.** This is a narrower claim than “intrinsic mesolimbic connectivity is age-invariant,” and the analyses below are designed to test it as such.

## 2. Methods

### 2.1 Participants and data

We used the NKI-Rockland sample (Enhanced Nathan Kline Institute Rockland Sample), a publicly available lifespan neuroimaging dataset. From the BIDS-organized release, we selected all resting-state fMRI scans with TR = 1400 ms (`task-rest_acq-1400_bold`). 428 such scans were available on the bucket for 307 unique participants; 417 of these yielded usable ROI time series and entered the analyses. Age, sex, and handedness were obtained from the accompanying `aws_links.csv` phenotypic file (the dataset's `participants.tsv` contains only sex/handedness; age is distributed separately). Participant age ranged from 9 to 83 years. We excluded no participants based on age; all scans with successful ROI time-series extraction were retained.

### 2.2 Regions of interest

The VTA and bilateral NAcc were defined as 4-mm-radius spheres in MNI space: VTA at (0, −16, −14) and NAcc at (±10, 12, −8). These coordinates were taken from Murty et al. (2014) to ensure comparability with the original study. The NAcc was defined bilaterally; the two hemispheres were averaged after correlation.

### 2.3 Functional coupling

For each scan, we extracted mean time series from the three spheres using `NiftiSpheresMasker` (nilearn), with detrending, z-standardization, and band-pass filtering (0.01–0.08 Hz; TR = 1.4 s). We computed the Pearson correlation between the VTA and each NAcc hemisphere using `ConnectivityMeasure`, averaged the two VTA–NAcc correlations, and applied the Fisher z-transformation. Scans with fewer than 30 volumes were excluded.

### 2.4 Statistical analysis

We tested age effects in three ways: (1) a simple linear Pearson correlation between age and Fisher-z coupling; (2) an inverse-age (1/age) model, replicating Murty et al. (2018), which assumes rapid early change followed by stabilization; and (3) a multiple regression controlling for sex. Because the NKI-Rockland sample includes multiple scans per participant, the primary analysis used **one observation per participant** (first scan, ordered by session label) to avoid pseudoreplication; we also report the full-scan analysis for completeness.

**Analysis selection and transparency.** This study was **not pre-registered**, and we therefore state the analysis plan explicitly rather than implying prior specification. The one-scan-per-participant rule was chosen for a methodological reason that does not depend on the results: with 81% of participants (249 of 307) contributing a single scan, treating repeated scans as independent would artificially inflate the effective degrees of freedom. The same rule is applied uniformly to every age model reported, including those that are not the primary analysis. We additionally report (i) the full-scan analysis, (ii) a linear mixed model with a random intercept per participant that uses all scans without discarding any (Section 3.7), and (iii) every sensitivity analysis we ran, in Figure 3. A reader who prefers a different primary specification can read the corresponding row of that figure; the conclusion does not depend on which row is chosen. Where we report a marginal effect (the full-scan inverse-age model, p = 0.068), we also report the sensitivity analyses that localise it.

### 2.5 Sensitivity analyses

To assess robustness, we repeated the analysis (i) excluding the youngest/oldest 5% of participants, (ii) excluding the youngest/oldest 10%, and (iii) excluding coupling outliers (|z| > 3 SD).

### 2.6 Head-motion quantification

The raw BIDS release used for the coupling estimates does not ship motion parameters. We therefore obtained framewise displacement (FD) from the C-PAC derivatives of the same dataset, which provide FD separately for each acquisition. Contrary to our initial expectation, FD **is** available for the 1400-ms resting scan from the C-PAC derivatives (see repository for the exact path).

Using this path, FD was retrieved for 307 participants. Of these, **250 had FD for their primary (first) scan**, and 304 of the 417 successfully processed scans could be matched to a session-level FD value. C-PAC did not process every session, so 57 participants have no FD for their primary scan; these participants are retained in the primary analysis but are necessarily excluded from the motion-controlled models. To confirm that this subset is not atypical, we compared included and excluded participants on age, sex, and coupling (Section 3.4).

### 2.7 Motion diagnostics

For a null result, the informative question is not only whether motion was controlled but whether it could have produced the null. We therefore report three diagnostics:

**D1**: the association between age and mean FD. If motion increases with age and reduces connectivity, an uncorrected analysis would be *biased toward* a negative age effect.

**D2**: the association between mean FD and coupling. If motion does not relate to coupling in this sample, then it cannot be mediating an age effect.

**D3**: comparison of included and excluded participants on age and coupling.

We then re-estimated the age effect in three ways: with mean FD as a covariate; after excluding participants above five FD thresholds (1.0, 0.7, 0.5, 0.3, and 0.2 mm); and in the full-scan dataset with mean FD as a covariate.

### 2.8 Equivalence testing and mixed-effects model

A non-significant test does not establish that an effect is absent, so we supplemented the null-hypothesis tests in two ways.

First, we performed **two one-sided tests (TOST)** for equivalence on the primary age–coupling correlation. TOST inverts the burden of proof: rather than asking whether the observed correlation differs from zero, it asks whether the data are compatible with the correlation lying inside a specified equivalence bound. We evaluated bounds from r = 0.05 to r = 0.20 and report the smallest bound that is supported at α = 0.05.

Second, we computed a **Bayes factor** for the point null (H0: ρ = 0) against a two-sided alternative, using the Fisher-z approximation with a Cauchy prior on the true z (scale 0.5 and 1.0; Wetzels & Wagenmakers, 2012). BF01 quantifies evidence *for* the absence of an effect, which is the quantity of interest for a replication of a null result.

Third, because repeated scans are available for only a minority of participants (58 of 307 overall, 22 of the 269 in the motion-matched subset), we fitted a **linear mixed model** with a random intercept per participant across all 304 scans with usable motion estimates (fixed effects: age, sex, mean FD). This uses the repeated scans without discarding them and without treating them as independent observations. The variance-ratio parameter was estimated by REML; the implementation was validated on simulated data with known intraclass correlations.

## 3. Results

### 3.1 Sample

417 scans from 307 unique participants were successfully processed (age 9–83 years). 249 participants contributed one scan, 26 contributed two, and 32 contributed three or more.

### 3.2 Age effects on VTA–NAcc resting-state coupling

Figure 1 shows the primary relationship; the fitted linear and inverse-age curves are both essentially flat across the 9–83 year range.

**Primary analysis (307 unique participants, one scan each).** There was no significant linear relationship between age and VTA–NAcc resting-state coupling: r = −0.042, p = 0.461. The inverse-age model (replicating Murty et al.) with sex as a covariate also showed no significant effect: β = −0.703, p = 0.339. Controlling for sex did not change the conclusion (linear model p = 0.814).

**Full-scan analysis (417 scans).** Results were consistent. Linear correlation: r = −0.001, p = 0.984. Age with sex as a covariate: β_age = 0.00040, p = 0.569, 95% CI = [−0.00096, 0.00176]. Inverse-age with sex: β = −1.135, p = 0.068 (marginal; without sex the same model gives β = −0.675, p = 0.258).

### 3.3 Sensitivity analyses

Excluding the youngest/oldest 5% of participants (n = 276, age 13–74) and 10% (n = 248, age 16–68) both eliminated any hint of an age effect: inverse-age p = 0.756 and p = 0.798, respectively; linear p = 0.405 and 0.474. Excluding coupling outliers (n = 305) likewise yielded no effect (linear p = 0.337; inverse-age p = 0.619).

The marginal inverse-age effect in the full-scan analysis (with sex; p = 0.068) was therefore driven by a small number of the youngest (<13 years) and oldest (>74 years) participants, combined with pseudoreplication from repeated scans; it was not a robust lifespan trend.

### 3.4 Head-motion diagnostics and control

Figure 2 presents the two motion diagnostics.

**D1 — motion increases with age.** Mean FD across the primary scans was 0.298 mm (SD = 0.157, range 0.045–1.016). Mean FD correlated positively with age: r = +0.177, p = 0.005, n = 250. Head motion is therefore a genuine age-related nuisance variable in this sample, and its typical effect (reducing apparent connectivity) would bias an uncorrected analysis toward a spurious negative age effect.

**D2 — motion is unrelated to coupling.** Mean FD was not associated with VTA–NAcc coupling: r = +0.056, p = 0.375, n = 250. The typical mechanism by which motion produces spurious age effects is not operating here: motion increases with age, but it does not suppress connectivity in this sample.

**D3 — included and excluded participants are comparable.** Participants with FD (n = 250) and without (n = 57) did not differ reliably in age (40.7 ± 19.7 vs 45.2 ± 21.6 years; Welch t = −1.44, p = 0.153) or in coupling (z = 0.17 ± 0.23 vs 0.08 ± 0.36; Welch t = 1.78, p = 0.079). Sex composition was similar (40.8% vs 36.8% male).

**Motion as a covariate.** Adding mean FD to the model left the age effect non-significant: linear age β = −0.00018, t = −0.225, p = 0.822; inverse-age β = −0.354, t = −0.470, p = 0.639. The mean FD coefficient itself was not significant in either model (p = 0.349 and p = 0.387).

**Excluding high-motion participants.** Figure 3 summarises the age effect across every model and sensitivity analysis; the partial correlation never exceeded |r| = 0.10. The age effect remained non-significant at every threshold tested:

| Mean FD threshold | n | linear p | age effect (controlling sex) p |
| --- | --- | --- | --- |
| ≤ 0.2 mm              | 68    | 0.884        | 0.417                              |
| ≤ 0.3 mm              | 148   | 0.611        | 0.324                              |
| ≤ 0.5 mm              | 231   | 0.925        | 0.603                              |
| ≤ 0.7 mm              | 242   | 0.750        | 0.874                              |
| ≤ 1.0 mm              | 249   | 0.656        | 0.975                              |

**Full-scan analysis with motion control.** Across 304 scans from 269 participants with mappable FD, the age effect was non-significant (linear r = 0.003, p = 0.955; with sex and mean FD, β_age = 0.00025, t = 0.325, p = 0.745).

### 3.5 Secondary finding: sex

Male participants showed higher VTA–NAcc resting-state coupling than female participants. In the full primary sample (n = 307) the effect was β = 0.073, t = 2.35, p = 0.019; in the motion-controlled subsample (n = 250) it was β = 0.065, t = 2.11, p = 0.036, and it remained significant with mean FD in the model (p = 0.036). This effect was not the target of the present study and is reported for completeness; we note it without a mechanistic interpretation.

![Figure 1](Fig1_age_vs_coupling.png)

**Figure 1.** Age and VTA–NAcc resting-state coupling, NKI-Rockland, n = 307.

![Figure 2](Fig2_motion_diagnostics.png)

**Figure 2.** Head-motion diagnostics. (A) Motion increases with age. (B) Motion is unrelated to coupling.

![Figure 3](Fig3_forest.png)

**Figure 3.** Partial correlation of age with coupling across all models and sensitivity analyses.

### 3.6 Equivalence testing

The 90% confidence interval for the primary age–coupling correlation (the interval used by TOST at α = 0.05) was [−0.136, +0.052]. The smallest equivalence bound supported by the data was **r = 0.136**: that is, we can reject age effects of |ρ| ≥ 0.14, but the present sample does not license excluding smaller effects (at a bound of r = 0.10, TOST p = 0.156).

The Bayes factor was BF01 = 8.49 with a Cauchy prior scale of 0.5 and BF01 = 16.74 with a scale of 1.0. By conventional criteria (Jeffreys, 1961) this is **moderate to strong evidence for the absence of an age effect**.

Taken together, the results are best described as: no detectable age effect, with medium-to-large effects excluded and small effects not excluded.

### 3.7 Linear mixed model

The random-intercept variance was estimated at approximately zero (variance ratio ≈ 0, implied ICC ≈ 0), so the mixed model reduced to ordinary least squares with participant-level degrees of freedom. The age effect remained non-significant (b = +0.00025, SE = 0.00075, t = 0.33, p = 0.745); the sex effect remained significant (b = +0.070, p = 0.020); mean FD was not significant (b = +0.111, p = 0.249). Simulation with matched design characteristics confirmed that the estimator recovers non-zero intraclass correlations (true ICC 0.20 → 0.21; 0.50 → 0.49), so the near-zero estimate is not an artefact of the number of repeated scans.

The near-zero ICC is itself informative: in this sample, repeated scans from the same participant contribute little clustering beyond the fixed effects, which means the pseudoreplication concern that motivates the one-scan-per-participant primary analysis is quantitatively small. It also suggests that scan-to-scan variability in VTA–NAcc coupling is substantial relative to stable between-person differences. Put differently, a single VTA–NAcc scan may capture a substantial share of measurement noise rather than a stable individual difference. Section 3.8 quantifies the reliability of the measure directly, and Section 4.6 sets out what this implies for the interpretation of the null.

![Figure 4](Fig4_robustness.png)

**Figure 4.** Robustness and reliability. (A) Agreement between the raw-BIDS and C-PAC preprocessed coupling estimates. (B) Test–retest reliability across BAS1–BAS2 session pairs. (C) Age–coupling correlation under twelve alternative ROI definitions.

### 3.8 Sensitivity analysis on fully preprocessed data

The analyses above used raw BIDS images, which carry no nuisance regression, CSF/white-matter regression, or volume scrubbing; FD entered only as a covariate. To test whether the conclusion depends on this, we recomputed coupling from the **C-PAC preprocessed derivatives** of the same acquisition, which are motion-corrected and nuisance-regressed, using the identical ROI definition and the identical filtering pipeline.

C-PAC did not process every session. Of the 409 scans planned, 171 sessions were absent from the derivatives (HTTP 404) and 220 were retrieved, covering 200 participants. Applying the same scan-selection rule as the primary analysis, first session per participant, yielded **n = 178** participants whose scan was available under both pipelines. This subset did not differ from the remaining 129 participants in age (42.0 ± 20.2 vs 40.9 ± 20.0 years; Welch t = 0.50, p = 0.620) or in sex composition (43.8% vs 34.9% male; z = 1.58, p = 0.115).

| analysis | n | r | p |
| --- | --- | --- | --- |
| 4-mm spheres, age only  | 178   | +0.045 | 0.553 |
| 4-mm spheres, age + sex | 178   | —      | 0.166 |
| + mean FD               | 178   | —      | 0.210 |

The age effect remains non-significant under full preprocessing. The point estimate is positive (r = +0.045) rather than negative as in the raw-BIDS analysis (r = −0.042), but the two are close in absolute magnitude and neither approaches significance.

Two further quantities help interpret the difference. First, the **correlation between the two coupling estimates is r = 0.281** (p = 0.0001, n = 178): two pipelines applied to the same scans agree only modestly. Second, the **test–retest reliability** of the measure can be estimated directly from participants scanned more than once. We report every estimate, because they differ substantially. Across **40 BAS1–BAS2 pairs** — the same session type as almost all primary scans — r = 0.279 and ICC(1,1) = 0.262. Across **all 182 pairs of sessions of any type**, which additionally includes influenza-challenge and treatment sessions, r = 0.111 and ICC(1,1) = 0.118. On the motion-matched subset used by the mixed model, the direct ICC is −0.04, consistent with that model's near-zero random-intercept variance.

The gap between the same-type and all-pairs figures reflects genuine heterogeneity across session types rather than estimation error: FLU sessions take place during an influenza challenge and TRT sessions during treatment, so including them inflates the apparent within-person variability of an otherwise stable measure. Because the primary analysis draws its single scan almost entirely from baseline sessions, **ICC ≈ 0.26 is the more relevant figure**; we nevertheless report the more conservative 0.12 as a lower bound.

The difference in sign between pipelines is therefore best read as noise rather than as a substantive contradiction: two measures with test–retest reliability in this range are expected to show limited agreement at the sample sizes available here. We do not claim that either point estimate is the correct one.

### 3.9 ROI definition sensitivity

To test whether the null depends on the spherical ROI definition, we recomputed coupling under twelve alternative definitions in the preprocessed data: sphere radii of 2, 3, 4, 5 and 6 mm; NAcc targets displaced to (±12, 12, −8); and VTA seeds displaced by ±2 mm along each of the three axes.

No variant produced a significant zero-order age–coupling correlation (all p > 0.26; point estimates ranged from r = +0.017 to +0.084, all positive). When sex was included as a covariate, the age term reached p < 0.05 in one of the twelve variants (VTA shifted −2 mm in x, p = 0.042). Given twelve comparisons, one nominal result of this magnitude is close to what would be expected by chance, and the corresponding zero-order correlation was not significant (r = +0.084, p = 0.267). We report the full set rather than only the primary ROI.

## 4. Discussion

We found that resting-state functional coupling between the VTA and NAcc shows **no detectable age effect** across the lifespan (9–83 years) in the NKI-Rockland sample. This replicates the resting-state null finding of Murty et al. (2018), extends it to a sample roughly 1.7× larger (307 vs. 180 participants) and an age range roughly 4× wider (74 vs. 20 years) that now includes middle and old age, and, for a null result, shows that the null is not an artefact of head motion.

### 4.1 Relation to prior developmental findings

Our result is convergent with, rather than independent of, two earlier developmental observations. Tomasi and Volkow (2014) did not find a significant VTA–NAcc difference between children/adolescents and adults, even though their abstract-level summary emphasises age-related increases in VTA–limbic connectivity more broadly; it is worth noting that in their data the *substantia nigra* seed did show significant developmental decreases with NAcc, which underlines the importance of distinguishing midbrain seeds. Morales et al. (2024) reported a VTA–NAcc age term of β = −0.006 (p = 0.83) across 4,025 longitudinal scans. The present study adds a formal replication framing, an extension into middle and old age, and a motion-controlled test.

### 4.2 The claim is edge-specific, not system-wide

Two findings constrain how broadly this result should be read.

First, within **the same NKI-Rockland data and the same 1400-ms acquisition**, O'Rawe and Leung (2022) found strong age-related degradation of caudate connectivity topography. Dorsal and ventral striatal territories therefore behave differently across the lifespan in this dataset. Our result should not be summarised as “striatal connectivity is age-invariant”; it concerns the VTA–NAcc edge specifically.

Second, **other VTA resting-state connections do change with age**: Peterson et al. (2017) reported age effects on VTA connectivity with temporal, parahippocampal, and cerebellar regions in adults. Stability is thus a property of this particular mesolimbic edge rather than of the VTA as a seed.

A plausible reading is that the ventral mesolimbic edge constitutes a relatively hard-wired motivational substrate whose *recruitment* is developmentally modulated, while dorsal striatal and other VTA connections undergo structural and functional reorganisation. Distinguishing these possibilities will require within-person longitudinal designs, which our cross-sectional analysis cannot provide.

### 4.3 Relation to the dopamine–ageing literature

Prestriatal dopamine transporter and D2 receptor availability decline with age, and Rieckmann et al. (2018) demonstrated that this decline co-occurs with caudate functional dedifferentiation, with both predicting memory change. Our finding of a stable VTA–NAcc edge sits alongside, rather than against, that framework: it suggests that molecular decline does not entail connectivity change uniformly across the striatum, and that ventral and dorsal territories may dissociate. We regard this dissociation as the most interesting implication of the present result, but we emphasise that it is an inference across studies rather than a within-study contrast, and we did not measure dopaminergic markers.

### 4.4 Head motion

Because this is a null result, motion control is central rather than incidental. Three observations jointly argue that motion did not produce the null. Motion increased with age (r = +0.18, p = 0.005), so any motion-induced suppression of connectivity would have *created* a negative age effect rather than masked one; yet motion was unrelated to coupling (r = +0.06, p = 0.38). The age effect remained non-significant with FD as a covariate (p = 0.82), after exclusion at five thresholds (all p > 0.32), and in the full-scan analysis (p = 0.75). Participants without FD did not differ reliably from those with it in age or coupling.

### 4.5 Specification choice and measurement reliability

A reader may notice that the primary specification (one scan per participant) yields the most straightforwardly null result, while the full-scan inverse-age model is marginal (p = 0.068). We think this concern deserves a direct answer rather than a footnote.

Three points are relevant. First, the one-scan-per-participant rule was applied to *every* model in the primary family, and it is not the specification that is most favourable to the null in the strict sense: the full-scan *linear* model is more null (p = 0.98) and the inverse-age model is less null. Second, the marginal full-scan effect does not survive the analyses designed to localise it. It disappears when the extreme age tails are trimmed (p = 0.76 and 0.80) and when mean FD is included (p = 0.75), which is what one expects if it is driven by a few extreme participants combined with pseudoreplication. Third, and most importantly, the mixed model in Section 3.7 uses **all** scans while treating the repeated-measures structure correctly, and it agrees with the primary analysis (age p = 0.75, ICC ≈ 0). A specification that uses every observation and requires no discarding of data reaches the same conclusion, which is stronger evidence than the primary analysis alone.

We nonetheless state plainly that had the full-scan inverse-age model been designated primary, the paper's headline claim would have had to be weaker. The numerical disagreement between specifications is small in absolute terms (partial r = −0.055 vs −0.014), and Figure 3 makes it visible rather than hiding it.

### 4.6 Limitations

**Motion was modelled, not removed.** The most important limitation is that the coupling estimates themselves were computed from raw BIDS images without nuisance regression, motion scrubbing, or CSF/white-matter regression; FD was used as a covariate and for exclusion. This addresses the age × motion confound but does not remove motion artefact from the connectivity estimates themselves. Section 3.8 reports the corresponding sensitivity analysis: coupling recomputed from fully preprocessed derivatives with confound regression remains non-significant (r = +0.045, p = 0.55), although only 178 of 307 participants could be included because C-PAC did not process every session.

**FD coverage.** FD was unavailable for the primary scan of 57 of 307 participants, because C-PAC did not process those sessions. Although these participants did not differ reliably from the analysed subset, the subsample is 81% of the full sample.

**Single acquisition.** We analysed only the 1400-ms TR resting-state scan. The NKI-Rockland sample also includes 645-ms multiband resting scans and 2500-ms scans; replicating across acquisitions would strengthen the result. This is especially relevant given that developmental resting-state fMRI is itself methodologically sensitive (Camacho et al., 2020).

**ROI definition.** Spherical ROIs are a crude approximation of the VTA and NAcc. Probabilistic atlases (e.g. Murty et al.'s VTA atlas; Harvard–Oxford subcortical atlas for NAcc) would provide more precise boundaries. We regard it as unlikely that finer ROI definition would reverse a null of this magnitude, but this should be verified.

**Cross-sectional inference.** The age effect here is cross-sectional; the longitudinal within-person trajectory is not estimated. Murty et al.'s original study was longitudinal; our cross-sectional replication is complementary but does not address within-person change.

**Null results and equivalence.** We report a non-significant age effect and, following TOST, do not claim equivalence at the level of small effects. With n = 307 and the observed effect size (r = −0.042), the 95% confidence interval for r is [−0.15, +0.07] and the smallest supported equivalence bound is r = 0.136. The data therefore exclude moderate and large associations but remain compatible with a small effect. The Bayes factor (BF01 = 8.5–16.7) provides positive evidence for the null, but positive evidence is not proof.

**Measurement reliability, and what the null can and cannot mean.** Single-scan VTA–NAcc coupling has, at best, fair reliability: ICC(1,1) = 0.262 across 40 same-type session pairs and 0.118 across all 182 pairs. Unreliable measurement attenuates observed correlations, so the true age effect may be larger than the observed r = −0.042 — by a factor of up to 1/√0.262 ≈ 2.0, or 1/√0.118 ≈ 2.9 on the more conservative estimate.

The consequence for the equivalence claim is substantial. Applying the same correction to the equivalence bound, the observed bound of r = 0.136 corresponds to a true bound of roughly 0.27 (same-type reliability) to 0.40 (all-pairs reliability). **The study therefore excludes large true age effects and provides only weak evidence against a medium one.** This is the single most important qualification of the present result, and it follows from the measurement properties of single-scan VTA–NAcc coupling rather than from the sample size. Future work should average across multiple scans per participant; at the observed reliability, raising reliability above 0.7 would require roughly six to ten scans per participant.

### 4.7 Conclusion

VTA–NAcc resting-state functional coupling shows **no detectable age effect** from childhood through old age, and this null survives head-motion control.

We use the phrase “no detectable age effect” rather than “stable” for two reasons. First, equivalence testing on the observed correlation excludes effects of |r| ≥ 0.14 but not smaller ones. Second, and more restrictively, the measure itself is only moderately reliable, so the same bound corresponds to a true effect of roughly |ρ| ≥ 0.27 after correction for attenuation. The appropriate claim is bounded twice over: the data are inconsistent with a large developmental effect on this edge, provide moderate-to-strong Bayesian evidence against any effect, and remain compatible with a medium effect. The 95% confidence interval for the primary correlation, [−0.15, +0.07], should be read alongside the conclusion rather than replaced by it.

Within those limits, the result supports the view that the tonic connectivity of this specific mesolimbic edge is developmentally conserved, and that age-related modulation of the mesolimbic system is expressed in motivated, task-evoked states and in dorsal striatal and other VTA connections rather than uniformly across the system.

## 5. Data and code availability

Data are from the NKI-Rockland sample, publicly available via the FCP-INDI S3 bucket (`s3://fcp-indi/data/Projects/RocklandSample/RawDataBIDSLatest`), the C-PAC derivatives (`s3://fcp-indi/data/Projects/RocklandSample/Outputs/C-PAC`), and `aws_links.csv` (phenotypic data). All analysis code, together with the derived subject-level tables needed to reproduce every statistic reported here without downloading imaging data, is available at **[repository DOI to be inserted]**. Analyses used Python 3.12 and nilearn; the statistics reported here were computed with a dependency-free implementation of the Student-t distribution (see repository), validated against reference values.

## References

- Camacho, M. C., Quiñones-Camacho, L. E., & Perlman, S. B. (2020). Does the child brain rest?: An examination and interpretation of resting cognition in developmental cognitive neuroscience. *NeuroImage, 212*, 116688. https://doi.org/10.1016/j.neuroimage.2020.116688
- Morales, A. M., Jones, S. A., Carlson, B., Kliamovich, D., Dehoney, J., Simpson, B. L., Dominguez-Savage, K. A., Hernandez, K. O., Lopez, D. A., Baker, F. C., Clark, D. B., Goldston, D. B., Luna, B., Nooner, K. B., Muller-Oehring, E. M., Tapert, S. F., Thompson, W. K., & Nagel, B. J. (2024). Associations between mesolimbic connectivity, and alcohol use from adolescence to adulthood. *Developmental Cognitive Neuroscience, 70*, 101478. https://doi.org/10.1016/j.dcn.2024.101478
- Murty, V. P., Shah, H., Montez, D., Foran, W., Calabro, F., & Luna, B. (2018). Age-related trajectories of functional coupling between the VTA and nucleus accumbens depend on motivational state. *Journal of Neuroscience, 38*(34), 7420–7427. https://doi.org/10.1523/JNEUROSCI.3508-17.2018
- Murty, V. P., Shermohammed, M., Smith, D. V., Carter, R. M., Huettel, S. A., & Adcock, R. A. (2014). Resting state networks distinguish human ventral tegmental area from substantia nigra. *NeuroImage, 100*, 580–589. https://doi.org/10.1016/j.neuroimage.2014.06.047
- Nooner, K. B., Colcombe, S. J., Tobe, R. H., Mennes, M., Benedict, M. M., Moreno, A. L., ... & Milham, M. P. (2012). The NKI-Rockland sample: A model for accelerating the pace of discovery science in psychiatry. *Frontiers in Neuroscience, 6*, 152. https://doi.org/10.3389/fnins.2012.00152
- O'Rawe, J. F., & Leung, H.-C. (2022). Topographic organization of the human caudate functional connectivity and age-related changes with resting-state fMRI. *Frontiers in Systems Neuroscience, 16*, 966433. https://doi.org/10.3389/fnsys.2022.966433
- Peterson, A. C., Zhang, S., Hu, S., Chao, H. H., & Li, C.-S. R. (2017). The effects of age, from young to middle adulthood, and gender on resting state functional connectivity of the dopaminergic midbrain. *Frontiers in Human Neuroscience, 11*, 52. https://doi.org/10.3389/fnhum.2017.00052
- Rieckmann, A., Johnson, K. A., Sperling, R. A., Buckner, R. L., & Hedden, T. (2018). Dedifferentiation of caudate functional connectivity and striatal dopamine transporter density predict memory change in normal aging. *Proceedings of the National Academy of Sciences, 115*(40), 10160–10165. https://doi.org/10.1073/pnas.1804641115
- Tomasi, D., & Volkow, N. D. (2014). Functional connectivity of substantia nigra and ventral tegmental area: Maturation during adolescence and effects of ADHD. *Cerebral Cortex, 24*(4), 935–944. https://doi.org/10.1093/cercor/bhs382
- Jeffreys, H. (1961). *Theory of probability* (3rd ed.). Oxford University Press.
- Wetzels, R., & Wagenmakers, E.-J. (2012). A default Bayesian hypothesis test for correlations and partial correlations. *Psychonomic Bulletin & Review, 19*(6), 1057–1064. https://doi.org/10.3758/s13423-012-0295-x
