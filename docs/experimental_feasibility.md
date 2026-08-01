# Experimental feasibility: predictions are not measurements

> **This package plans hypothesis tests. A computed spectrum or Fourier peak is
> not evidence that a biological Davydov soliton exists.**

## Evidence layers and validation path

1. **Hypothesis:** a coupled excitation may localize and affect folding.
2. **Theoretical prediction:** a specified model predicts dynamics; unmeasured
   parameters remain assumptions.
3. **Simulated observable:** a response model predicts a 2D-IR feature.
4. **Measured quantity:** an instrument records signal, noise, isotope shift,
   concentration dependence, and control contrast.

```text
Stage A: labeled synthetic peptide
  │ measured sensitivity + isotope/control criteria pass
  ▼
Stage B: purified stalled RNCs
  │ reporter concentration exceeds the empirical threshold
  ▼
Stage C: active translation
    synchronization + occupancy + sensitivity thresholds all pass
```

Stage A uses isotope-edited backbone labels, confined/unconfined and
sequence-scrambled controls, stiffness-changing mutations, temperature series,
and measured concentration-dependent limits. Stage B adds empty ribosomes,
nonfolding and released-chain controls, multiple chain lengths, and independent
biological purifications. Stage C adds no-mRNA, initiation-only, defined-stall,
puromycin-release, and predicted sequence-perturbation controls.

## Population budget (SI units)

Inputs are RNC concentration (mol/L), launch rate (s⁻¹), lifetime (s),
synchronization and competence fractions (0–1), reporter count, labeling
efficiency (0–1), and empirical detection threshold (mol/L). There is no silent
unit conversion.

`f_active = min(1, k_launch * tau_state)`

`C_active = C_RNC * f_active * f_synchronized * f_competent`

Reporter concentration additionally multiplies reporter count and labeling
efficiency. Heterodyne signal is assumed approximately **linear** in molecular
concentration; no quadratic concentration enhancement is used. Shortfall,
required concentration, and required lifetime expose zero-occupancy cases as
`Inf` rather than dividing by zero.

### Worked no-go: active translation

The editable default uses 100 nM (`100e-9` mol/L), 10 s⁻¹ launch, 1 ps
(`1e-12` s) lifetime, 0.1 synchronization, 0.5 competence, two reporters, 0.8
labeling, and a 1 μM (`1e-6` mol/L) threshold. Its active fraction is `1e-11`
and effective reporter concentration `8e-19` mol/L: a no-go for ordinary
ensemble averaging.

### Worked feasible labeled peptide

The synthetic-peptide default uses 500 μM (`5e-4` mol/L), two reporters, 0.9
labeling, and a measured 10 μM (`1e-5` mol/L) threshold. Its effective reporter
concentration is 900 μM. This planning example is instrument-feasible, while
solubility, aggregation, path length, and response still require measurement.

## Separate clocks

`CoherenceTime`, `WaitingTime`, `DetectionTime`, and
`BiologicalReactionTime` accept seconds and reject empty, negative, non-finite,
repeated, or non-monotonic arrays. A femtosecond 2D-IR waiting-time scan does
**not** synchronize peptide-bond formation after a biological trigger.

## Uncertainty, power, and beat analysis

`UniformRange` with `monte_carlo_feasibility` gives seeded sampling, a median,
95% interval, go/threshold probabilities, and sensitivity ranking.
`sensitivity_analysis` returns a plotting-ready grid over concentration,
lifetime, synchronization, reporters, labeling, averages, and empirical noise.
Its √N scaling assumes independent averages and must be replaced when measured
instrument behavior differs. It makes no acquisition-time claim.

`analyze_beats` validates uniform sampling and Nyquist, optionally detrends and
windows, reports resolution, warns about short traces, and compares a sinusoid
with a no-oscillation model using ΔBIC. `candidate_frequency` and
`spectral_evidence` are neutral outputs: a Fourier peak alone does not establish
quantum coherence or a soliton.

## Decisions and limitations

`recommend_stage` returns criteria for concentration/noise, isotope resolution,
control effect size, attainable concentration, and plausible lifetime. A passing
feature is *consistent with the predicted coupled mode* and may support
progression; it is insufficient to distinguish the mechanism from alternatives.
The package does not model heating, aggregation, spectral congestion, detector
nonlinearity, correlated noise, multiple testing, or a validated null ensemble.
