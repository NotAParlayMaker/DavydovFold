"""Parameters used to test whether a transient state is detectable in bulk 2D-IR.

All concentrations are molar and all rates/times are SI. `detection_limit_M`
should be replaced by an empirically measured limit for the planned instrument,
reporter, path length, and acquisition protocol.
"""
Base.@kwdef struct FeasibilityParams
    rnc_concentration_M::Float64 = 100e-9
    elongation_rate_s::Float64 = 20.0
    state_lifetime_s::Float64 = 200e-12
    detection_limit_M::Float64 = 1e-6
    synchronized_fraction::Float64 = 1.0
    labels_per_rnc::Int = 1
end

function _validate_feasibility(p::FeasibilityParams)
    values = (p.rnc_concentration_M, p.elongation_rate_s, p.state_lifetime_s,
              p.detection_limit_M, p.synchronized_fraction)
    all(isfinite, values) || throw(ArgumentError("feasibility parameters must be finite"))
    p.rnc_concentration_M >= 0 || throw(ArgumentError("rnc_concentration_M must be nonnegative"))
    p.elongation_rate_s >= 0 || throw(ArgumentError("elongation_rate_s must be nonnegative"))
    p.state_lifetime_s >= 0 || throw(ArgumentError("state_lifetime_s must be nonnegative"))
    p.detection_limit_M > 0 || throw(ArgumentError("detection_limit_M must be positive"))
    0 <= p.synchronized_fraction <= 1 ||
        throw(ArgumentError("synchronized_fraction must be between zero and one"))
    p.labels_per_rnc >= 1 || throw(ArgumentError("labels_per_rnc must be at least one"))
    p
end

"""Estimate the steady-state population of a short-lived, elongation-triggered state.

The occupancy is modeled as the probability that at least one elongation event
occurred during the state's lifetime, `1 - exp(-rate*lifetime)`. This makes the
usually hidden population-budget assumption explicit and avoids treating every
ribosome as continuously occupied by the transient state.
"""
function estimate_transient_population(p::FeasibilityParams=FeasibilityParams())
    _validate_feasibility(p)
    occupancy = -expm1(-p.elongation_rate_s * p.state_lifetime_s)
    active = p.rnc_concentration_M * occupancy * p.synchronized_fraction
    reporter = active * p.labels_per_rnc
    margin = reporter / p.detection_limit_M

    denominator = occupancy * p.synchronized_fraction * p.labels_per_rnc
    required_rnc = denominator > 0 ? p.detection_limit_M / denominator : Inf

    available_reporter = p.rnc_concentration_M * p.synchronized_fraction * p.labels_per_rnc
    target_occupancy = available_reporter > 0 ? p.detection_limit_M / available_reporter : Inf
    required_lifetime = if target_occupancy <= 0
        0.0
    elseif target_occupancy >= 1 || p.elongation_rate_s == 0
        Inf
    else
        -log1p(-target_occupancy) / p.elongation_rate_s
    end

    (occupancy_fraction=occupancy,
     active_state_concentration_M=active,
     reporter_concentration_M=reporter,
     detection_margin=margin,
     feasible=margin >= 1,
     required_rnc_concentration_M=required_rnc,
     required_state_lifetime_s=required_lifetime)
end

"""Return a go/no-go assessment and the least ambitious defensible test stage."""
function assess_2dir_feasibility(p::FeasibilityParams=FeasibilityParams())
    population = estimate_transient_population(p)
    recommendations = String[]

    if !population.feasible
        push!(recommendations,
              "Do not begin with active translation; the predicted reporter population is below the stated detection limit.")
        push!(recommendations,
              "Measure the instrument-specific detection limit with the chosen isotope or non-natural vibrational reporter.")
    end
    p.synchronized_fraction < 0.5 &&
        push!(recommendations, "Improve biochemical synchronization or model the unsynchronized population explicitly.")
    p.labels_per_rnc == 1 &&
        push!(recommendations, "Evaluate multiple site-specific reporters before relying on unlabeled amide-I contrast.")

    recommended_stage = if population.feasible
        :active_translation
    elseif p.rnc_concentration_M >= p.detection_limit_M
        :purified_stalled_rnc
    else
        :synthetic_peptide
    end

    (; population...,
       status=population.feasible ? :go : :no_go,
       recommended_stage=recommended_stage,
       recommendations=recommendations)
end

"""Return a staged validation ladder with explicit experimental gates."""
function staged_validation_plan()
    [
        (stage=:synthetic_peptide,
         system="20-40 residue peptide with site-specific vibrational labels",
         primary_endpoint="A pre-registered isotope-dependent band and coupling pattern",
         controls=["scrambled sequence", "stiffness-changing mutant", "unconfined peptide"],
         gate="Observed frequency, linewidth, and perturbation response agree with forward simulations"),
        (stage=:purified_stalled_rnc,
         system="Homogeneous stalled ribosome-nascent-chain complexes at defined chain lengths",
         primary_endpoint="A chain-length-dependent reporter feature above empty-ribosome background",
         controls=["empty ribosome", "nonfolding nascent chain", "released nascent chain"],
         gate="Signal reproduces across preparations and follows the predicted label and chain-length dependence"),
        (stage=:active_translation,
         system="Interleaved translating and matched nontranslating reactions",
         primary_endpoint="A translation-dependent differential spectrum that scales with active RNC concentration",
         controls=["no mRNA", "defined stall", "puromycin release", "model-suppressing sequence variant"],
         gate="Blinded model comparison rejects thermal, population-transfer, and ordinary structural alternatives")
    ]
end

"""Estimate a relative heterodyne signal using linear concentration scaling.

This is a transparent toy model, not an instrument calibration. Concentration
is in molar units. `active_fraction` can be supplied from
`estimate_transient_population`; the function deliberately does not apply an
N-squared molecular-concentration enhancement.
"""
function predict_signal_amplitude(rnc_concentration::Real, pulse_energy::Real;
                                  efficiency::Real=0.15, noise_floor::Real=1e-6,
                                  active_fraction::Real=1.0, labels_per_rnc::Integer=1)
    rnc_concentration >= 0 && pulse_energy >= 0 || throw(ArgumentError("inputs must be nonnegative"))
    efficiency >= 0 && noise_floor > 0 || throw(ArgumentError("efficiency must be nonnegative and noise_floor positive"))
    0 <= active_fraction <= 1 || throw(ArgumentError("active_fraction must be between zero and one"))
    labels_per_rnc >= 1 || throw(ArgumentError("labels_per_rnc must be at least one"))
    reporter_concentration = rnc_concentration * active_fraction * labels_per_rnc
    signal = efficiency * reporter_concentration * pulse_energy^1.5
    (reporter_concentration=reporter_concentration, signal=signal, snr=signal/noise_floor)
end

"""Design uniformly spaced 2D-IR waiting times and a four-step phase cycle."""
function design_pulse_sequence(desired_t2_range; step::Real=0.05)
    isfinite(step) && step > 0 || throw(ArgumentError("step must be finite and positive"))
    lo,hi=extrema(desired_t2_range); waits=collect(lo:step:hi)
    (waiting_times=waits, phase_cycle=[0,pi/2,pi,3pi/2], coherence_delays=(-0.5,0.5))
end

"""Design a separate biological reaction-time scan.

Reaction time is intentionally separated from the 2D-IR waiting time so that
millisecond-to-second translation synchronization is not confused with the
femtosecond-to-picosecond pulse sequence.
"""
function design_reaction_scan(desired_reaction_range; step::Real)
    isfinite(step) && step > 0 || throw(ArgumentError("step must be finite and positive"))
    lo, hi = extrema(desired_reaction_range)
    isfinite(lo) && isfinite(hi) || throw(ArgumentError("reaction-time bounds must be finite"))
    lo >= 0 || throw(ArgumentError("reaction times must be nonnegative"))
    (reaction_times=collect(lo:step:hi), units=:seconds)
end

"""Score a predicted beat in a two-column, uniformly sampled CSV trace."""
function analyze_beats(path::AbstractString; predicted_beat_freq::Real=1.2)
    isfinite(predicted_beat_freq) && predicted_beat_freq >= 0 ||
        throw(ArgumentError("predicted_beat_freq must be finite and nonnegative"))
    data=readdlm(path, ',', Float64)
    size(data, 2) >= 2 || throw(ArgumentError("trace must contain at least two columns"))
    t=data[:,1]; y=data[:,2].-mean(data[:,2])
    length(t) >= 4 || throw(ArgumentError("at least four samples are required"))
    dt = sampling_interval(t; name="sample times")
    all(isfinite, y) || throw(ArgumentError("signal values must contain only finite values"))
    f=FFTW.rfftfreq(length(y),inv(dt)); power=abs2.(FFTW.rfft(y))
    i=argmin(abs.(f.-predicted_beat_freq)); baseline=median(power[2:end])
    (frequency=f[i], likelihood_ratio=power[i]/max(baseline,eps()), power=power[i])
end

const CODONS=Dict('A'=>["GCT","GCC","GCA","GCG"], 'G'=>["GGT","GGC","GGA","GGG"],
 'L'=>["TTA","TTG","CTT","CTC","CTA","CTG"], 'P'=>["CCT","CCC","CCA","CCG"],
 'F'=>["TTT","TTC"], 'M'=>["ATG"], 'K'=>["AAA","AAG"], 'E'=>["GAA","GAG"],
 'D'=>["GAT","GAC"], 'V'=>["GTT","GTC","GTA","GTG"], 'I'=>["ATT","ATC","ATA"],
 'R'=>["CGT","CGC","CGA","CGG","AGA","AGG"], 'H'=>["CAT","CAC"],
 'N'=>["AAT","AAC"], 'Q'=>["CAA","CAG"], 'S'=>["TCT","TCC","TCA","TCG","AGT","AGC"],
 'T'=>["ACT","ACC","ACA","ACG"], 'C'=>["TGT","TGC"], 'W'=>["TGG"], 'Y'=>["TAT","TAC"])

"""Select synonymous codons using the sign of a desired local translation profile.

Positive entries select the first (nominally fast) codon and non-positive entries
the last (nominally slow) codon. Invalid residue codes are rejected.
"""
function codon_optimizer(sequence::AbstractString,target_profile::AbstractVector{<:Real})
    length(sequence)==length(target_profile) || throw(DimensionMismatch("profile must match sequence"))
    [begin
         choices=get(CODONS,uppercase(aa), nothing)
         isnothing(choices) && throw(ArgumentError("unsupported amino acid: $aa"))
         choices[target_profile[i]>0 ? 1 : length(choices)]
     end for (i,aa) in pairs(sequence)]
end
