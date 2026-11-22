# Staged Experimental Validation Program

This document converts the proposed 2D-IR test into a sequence of falsifiable,
resource-aware experiments. It does not assume that a simulated spectral feature
is experimentally observable merely because it exists in the model.

## 1. Population budget before spectroscopy

For a transient launched by elongation, the occupied fraction is modeled as

```text
occupancy = 1 - exp(-elongation_rate × state_lifetime)
```

and the reporter concentration is

```text
reporter concentration = RNC concentration × occupancy × synchronized fraction × labels per RNC
```

Run `estimate_transient_population` before choosing a sample format. The direct
active-translation experiment is a no-go whenever the calculated reporter
concentration is below the empirically measured detection limit.

The software reports both the required RNC concentration and the required state
lifetime. Infinite required lifetime means the chosen total reporter pool cannot
reach the detection threshold even at complete occupancy.

## 2. Three-stage validation ladder

### Stage A: labeled synthetic peptide

**System:** A 20-40 residue peptide containing one or more site-specific
vibrational labels.

**Primary endpoint:** A pre-registered isotope-dependent band, cross-peak, or
beat component whose frequency and perturbation response agree with a forward
spectral simulation.

**Minimum controls:**

- scrambled sequence
- stiffness-changing mutation
- labeled versus unlabeled peptide
- confined versus unconfined peptide
- pulse-energy and temperature series

**Go criterion:** The candidate component is reproducible across preparations,
moves as predicted under isotope substitution, and beats a prespecified
alternative model.

### Stage B: purified stalled RNCs

**System:** Homogeneous stalled ribosome-nascent-chain complexes at several
defined chain lengths.

**Primary endpoint:** A chain-length-dependent reporter feature that exceeds
empty-ribosome and released-chain backgrounds.

**Minimum controls:**

- empty ribosome
- nonfolding nascent chain
- released nascent chain
- model-suppressing sequence variant
- independent RNC preparations

**Go criterion:** The signal follows the predicted label position and chain
length, survives batch replication, and cannot be explained by sample
composition or ordinary secondary-structure formation.

### Stage C: active translation

**System:** Interleaved translating and matched nontranslating reactions.

**Primary endpoint:** A differential signal that scales with active RNC
concentration and follows a model-specific sequence perturbation.

**Minimum controls:**

- no-mRNA reaction
- defined biochemical or genetic stall
- puromycin release
- translation-rate and RNC-occupancy measurements
- model-suppressing sequence variant

**Go criterion:** Blinded model comparison rejects thermal response, ordinary
vibrational coherence, population transfer, and bulk structural changes as
sufficient explanations.

## 3. Keep the clocks separate

A 2D-IR experiment contains an internal coherence time and waiting time. A
translation experiment also has a biological reaction time. These axes must not
be conflated.

Use `design_pulse_sequence` for the internal spectroscopic waiting-time axis and
`design_reaction_scan` for the biological trigger-to-probe axis. Any claim of a
sub-picosecond biological rise requires a trigger that prepares the molecular
event with comparable synchronization.

## 4. Pre-registration fields

Before collecting the decisive dataset, record:

- sample concentration and active RNC fraction
- label identity, number, and position
- measured noise-equivalent concentration or absorbance
- path length, pulse energy, repetition rate, and acquisition time
- candidate peak windows with uncertainty bounds
- expected sign and phase of cross-peaks
- primary endpoint and exclusion criteria
- null, classical, and coupled-mode comparison models
- number of independent biological and technical replicates
- analysis code version and random seed

Exact peak coordinates should come from an explicit forward model and should be
reported with uncertainty ranges rather than as universal constants.

## 5. Interpretation guardrails

A long-lived oscillation is not by itself evidence of a soliton or topology.
Report a feature as consistent with the proposed coupled state only after
ordinary vibrational coherence, heating, spectral diffusion, and population
transfer have been quantitatively compared.

Claims of topological protection require a defined invariant and a demonstrated
robustness property. Claims of unambiguous detection require the prespecified
alternative models to fail, not merely the preferred model to fit.

## 6. Recommended first milestone

The first credible milestone is:

> A pre-specified, isotope-dependent vibrational feature appears in a confined
> model peptide and in purified RNCs of the appropriate length; its frequency,
> coupling pattern, temperature dependence, and mutation response agree with
> blinded forward simulations.

Passing that milestone justifies the active-translation experiment. Failing it
provides a clean falsification or a concrete reason to revise the model.
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
