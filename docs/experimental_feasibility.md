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
