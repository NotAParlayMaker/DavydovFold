"""Estimate signal and heterodyne SNR from concentration and pulse energy."""
function predict_signal_amplitude(rnc_concentration::Real, pulse_energy::Real;
                                  efficiency::Real=0.15, noise_floor::Real=1e-6)
    rnc_concentration >= 0 && pulse_energy >= 0 || throw(ArgumentError("inputs must be nonnegative"))
    efficiency >= 0 && noise_floor > 0 || throw(ArgumentError("efficiency must be nonnegative and noise_floor positive"))
    signal=efficiency*rnc_concentration*pulse_energy^1.5
    (signal=signal, snr=signal/noise_floor)
end

"""Create waiting times and a four-step phase cycle for a requested range."""
function design_pulse_sequence(desired_t2_range; step::Real=0.05)
    step > 0 || throw(ArgumentError("step must be positive"))
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
"""A validated time axis. `values_s` are seconds and must be finite, nonnegative,
and strictly increasing. The type parameter distinguishes physically different clocks."""
struct ExperimentTimeAxis{K}
    values_s::Vector{Float64}
    function ExperimentTimeAxis{K}(values_s::AbstractVector{<:Real}) where K
        v = collect(Float64, values_s)
        isempty(v) && throw(ArgumentError("time axis must not be empty"))
        all(isfinite, v) && all(>=(0), v) || throw(ArgumentError("times must be finite and nonnegative seconds"))
        all(diff(v) .> 0) || throw(ArgumentError("times must be strictly increasing"))
        new{K}(v)
    end
end

"""2D-IR coherence-time axis in seconds."""
CoherenceTime(values_s) = ExperimentTimeAxis{:coherence}(values_s)
"""2D-IR population/waiting-time axis in seconds; not a biological synchronization clock."""
WaitingTime(values_s) = ExperimentTimeAxis{:waiting}(values_s)
"""Signal-detection-time axis in seconds."""
DetectionTime(values_s) = ExperimentTimeAxis{:detection}(values_s)
"""Biological reaction time after a trigger or stall release, in seconds."""
BiologicalReactionTime(values_s) = ExperimentTimeAxis{:biological}(values_s)

"""Editable inputs for a transient-state population budget. Concentrations are mol/L."""
Base.@kwdef struct PopulationBudget
    rnc_concentration_mol_L::Float64
    launch_rate_s_inv::Float64
    state_lifetime_s::Float64
    synchronized_fraction::Float64 = 1.0
    competent_fraction::Float64 = 1.0
    reporters_per_active_state::Float64 = 1.0
    labeling_efficiency::Float64 = 1.0
    detection_threshold_mol_L::Float64
    signal_model::Symbol = :linear
end

function _validate(p::PopulationBudget)
    for (name, x) in ((:rnc_concentration_mol_L,p.rnc_concentration_mol_L),
                      (:launch_rate_s_inv,p.launch_rate_s_inv), (:state_lifetime_s,p.state_lifetime_s),
                      (:reporters_per_active_state,p.reporters_per_active_state),
                      (:detection_threshold_mol_L,p.detection_threshold_mol_L))
        isfinite(x) && x >= 0 || throw(ArgumentError("$name must be finite and nonnegative"))
    end
    p.detection_threshold_mol_L > 0 || throw(ArgumentError("detection threshold must be positive"))
    for (name, x) in ((:synchronized_fraction,p.synchronized_fraction),
                      (:competent_fraction,p.competent_fraction), (:labeling_efficiency,p.labeling_efficiency))
        isfinite(x) && 0 <= x <= 1 || throw(ArgumentError("$name must be in [0, 1]"))
    end
    p.reporters_per_active_state > 0 || throw(ArgumentError("reporter count must be positive"))
    p.signal_model == :linear || throw(ArgumentError("only the documented :linear heterodyne concentration model is supported"))
end

"""Calculate instantaneous occupancy and reporter concentration using a linear
heterodyne concentration model. No simulated feature is interpreted as evidence."""
function population_budget(p::PopulationBudget)
    _validate(p)
    active_fraction = min(1.0, p.launch_rate_s_inv * p.state_lifetime_s)
    availability = p.synchronized_fraction * p.competent_fraction
    active = p.rnc_concentration_mol_L * active_fraction * availability
    reporter = active * p.reporters_per_active_state * p.labeling_efficiency
    shortfall = max(0.0, p.detection_threshold_mol_L - reporter)
    scale = active_fraction * availability * p.reporters_per_active_state * p.labeling_efficiency
    required_rnc = scale == 0 ? Inf : p.detection_threshold_mol_L / scale
    lifetime_scale = p.rnc_concentration_mol_L * p.launch_rate_s_inv * availability *
                     p.reporters_per_active_state * p.labeling_efficiency
    raw_required_lifetime = lifetime_scale == 0 ? Inf : p.detection_threshold_mol_L / lifetime_scale
    required_lifetime = raw_required_lifetime <= inv(max(p.launch_rate_s_inv, eps())) ? raw_required_lifetime : Inf
    feasible = reporter >= p.detection_threshold_mol_L
    explanation = feasible ?
        "Effective reporter concentration meets the empirical threshold; this supports measurement, not mechanistic proof." :
        "Effective reporter concentration is below the empirical threshold; test a more concentrated or longer-lived system first."
    (; active_fraction, active_state_concentration_mol_L=active,
       effective_reporter_concentration_mol_L=reporter, concentration_shortfall_mol_L=shortfall,
       required_rnc_concentration_mol_L=required_rnc, required_state_lifetime_s=required_lifetime,
       feasibility_status=feasible ? :go : :no_go, explanation)
end

"""Return an editable named feasibility scenario; defaults are planning assumptions, not constants."""
function feasibility_scenario(name::Symbol)
    name == :synthetic_peptide && return PopulationBudget(rnc_concentration_mol_L=5e-4,
        launch_rate_s_inv=1.0, state_lifetime_s=1.0, detection_threshold_mol_L=1e-5,
        reporters_per_active_state=2.0, labeling_efficiency=0.9)
    name == :stalled_rnc && return PopulationBudget(rnc_concentration_mol_L=5e-6,
        launch_rate_s_inv=1.0, state_lifetime_s=1.0, synchronized_fraction=0.8,
        competent_fraction=0.5, reporters_per_active_state=2.0, labeling_efficiency=0.8,
        detection_threshold_mol_L=1e-6)
    name == :active_translation && return PopulationBudget(rnc_concentration_mol_L=100e-9,
        launch_rate_s_inv=10.0, state_lifetime_s=1e-12, synchronized_fraction=0.1,
        competent_fraction=0.5, reporters_per_active_state=2.0, labeling_efficiency=0.8,
        detection_threshold_mol_L=1e-6)
    throw(ArgumentError("unknown scenario $name; use :synthetic_peptide, :stalled_rnc, or :active_translation"))
end

const STAGE_CONTROLS = Dict(
 :A => ["isotope-edited backbone labels", "confined versus unconfined samples", "sequence-scrambled controls", "stiffness-changing mutations", "temperature series", "measured concentration-dependent detection limits"],
 :B => ["empty ribosome control", "nonfolding nascent-chain control", "multiple nascent-chain lengths", "released-chain control", "independently purified biological replicates"],
 :C => ["translating sample", "no-mRNA control", "initiation-only control", "defined stall control", "puromycin-release control", "sequence perturbation predicted to alter the signal"])

"""Recommend the earliest defensible validation stage and return machine-readable criteria."""
function recommend_stage(result; isotope_shift_resolvable::Bool=false, control_effect_size::Real=0,
                         minimum_effect_size::Real=1, attainable_concentration_mol_L::Real=1e-3,
                         plausible_lifetime_s::Tuple{<:Real,<:Real}=(0.0, Inf),
                         synchronization_threshold::Real=0.5, occupancy_threshold::Real=0.1,
                         inputs::Union{Nothing,PopulationBudget}=nothing)
    minimum_effect_size >= 0 && attainable_concentration_mol_L >= 0 || throw(ArgumentError("criteria bounds must be nonnegative"))
    threshold_met = result.feasibility_status == :go
    attainable = result.required_rnc_concentration_mol_L <= attainable_concentration_mol_L
    lifetime_plausible = plausible_lifetime_s[1] <= result.required_state_lifetime_s <= plausible_lifetime_s[2]
    contrast = control_effect_size >= minimum_effect_size
    criteria = (; minimum_reporter_concentration_reached=threshold_met,
        predicted_feature_exceeds_empirical_noise_threshold=threshold_met,
        isotope_shift_resolvable, control_contrast_sufficient=contrast,
        required_sample_concentration_attainable=attainable, required_lifetime_plausible=lifetime_plausible)
    stage = :A
    if inputs !== nothing && threshold_met
        occupancy = min(1.0, inputs.launch_rate_s_inv * inputs.state_lifetime_s)
        stage = inputs.synchronized_fraction >= synchronization_threshold && occupancy >= occupancy_threshold &&
                inputs.competent_fraction >= occupancy_threshold ? :C : :B
    elseif threshold_met
        stage = :B
    end
    (; stage, criteria, controls=STAGE_CONTROLS[stage],
       interpretation="A passing result is consistent with the predicted coupled mode and supports progression; it is insufficient alone to distinguish the proposed mechanism from alternatives.")
end

"""Uniform uncertainty interval sampled without additional dependencies."""
struct UniformRange
    low::Float64
    high::Float64
    function UniformRange(low::Real, high::Real)
        isfinite(low) && isfinite(high) && 0 <= low <= high || throw(ArgumentError("range must be finite, nonnegative, and ordered"))
        new(low, high)
    end
end
_draw(rng, x::Real) = Float64(x)
_draw(rng, x::UniformRange) = x.low + rand(rng) * (x.high - x.low)

"""Monte Carlo population feasibility with a reproducible seed and Spearman-like
absolute rank-correlation sensitivity ranking. Uncertain fields accept `UniformRange`."""
function monte_carlo_feasibility(; samples::Integer=10_000, seed::Integer=42, kwargs...)
    samples > 1 || throw(ArgumentError("samples must exceed one"))
    rng=MersenneTwister(seed); names=collect(keys(kwargs)); draws=Dict(n => Float64[] for n in names)
    active=Float64[]; decisions=Bool[]; threshold=Float64[]
    for _ in 1:samples
        vals=Dict(n => _draw(rng, kwargs[n]) for n in names)
        foreach(n -> push!(draws[n], vals[n]), names)
        p=PopulationBudget(; (n => vals[n] for n in names)...)
        r=population_budget(p); push!(active,r.active_state_concentration_mol_L)
        push!(decisions,r.feasibility_status == :go); push!(threshold,p.detection_threshold_mol_L)
    end
    sorted=sort(active); q(x)=sorted[clamp(round(Int, x*(samples-1))+1,1,samples)]
    rank(v)=sortperm(sortperm(v)); ar=rank(active)
    sensitivity=sort([(parameter=n, score=length(unique(draws[n])) < 2 ? 0.0 :
        abs(cor(Float64.(rank(draws[n])),Float64.(ar)))) for n in names], by=x->-x.score)
    reporters=[population_budget(PopulationBudget(; (n => draws[n][i] for n in names)...)).effective_reporter_concentration_mol_L for i in 1:samples]
    (; median_active_concentration_mol_L=median(active), credible_interval_mol_L=(q(0.025),q(0.975)),
       probability_exceeding_detection_threshold=mean(reporters .>= threshold),
       probability_go=mean(decisions), sensitivity_ranking=sensitivity, seed, samples)
end

"""Evaluate a Cartesian sensitivity grid and return plot-ready named tuples."""
function sensitivity_analysis(base::PopulationBudget; rnc_concentrations_mol_L=[base.rnc_concentration_mol_L],
 lifetime_s=[base.state_lifetime_s], synchronization=[base.synchronized_fraction],
 reporter_counts=[base.reporters_per_active_state], labeling_efficiencies=[base.labeling_efficiency],
 number_of_averages=[1], empirical_noise_floors_mol_L=[base.detection_threshold_mol_L])
    rows=NamedTuple[]
    for c in rnc_concentrations_mol_L, life in lifetime_s, sync in synchronization, reporters in reporter_counts,
        label in labeling_efficiencies, navg in number_of_averages, noise in empirical_noise_floors_mol_L
        navg isa Integer && navg > 0 || throw(ArgumentError("number of averages must be a positive integer"))
        # Independent-noise assumption; replace with measured instrument scaling when available.
        threshold=Float64(noise)/sqrt(navg)
        p=PopulationBudget(c,base.launch_rate_s_inv,life,sync,base.competent_fraction,reporters,label,threshold,base.signal_model)
        push!(rows,(; rnc_concentration_mol_L=Float64(c), state_lifetime_s=Float64(life),
            synchronization_fraction=Float64(sync), reporters_per_active_state=Float64(reporters),
            labeling_efficiency=Float64(label), number_of_averages=navg,
            empirical_noise_floor_mol_L=Float64(noise), result=population_budget(p)))
    end
    rows
end

"""Estimate heterodyne SNR using an explicitly linear concentration model."""
function predict_signal_amplitude(concentration_mol_L::Real, pulse_energy::Real;
                                  efficiency::Real=0.15, noise_floor::Real=1e-6)
    concentration_mol_L >= 0 && pulse_energy >= 0 || throw(ArgumentError("inputs must be nonnegative"))
    efficiency >= 0 && noise_floor > 0 || throw(ArgumentError("efficiency must be nonnegative and noise_floor positive"))
    signal=efficiency*concentration_mol_L*pulse_energy
    (signal=signal, snr=signal/noise_floor, concentration_model=:linear)
end

"""Design uniformly spaced 2D-IR waiting times (in the caller's stated unit)."""
function design_pulse_sequence(desired_t2_range; step::Real=0.05)
    isfinite(step) && step > 0 || throw(ArgumentError("step must be finite and positive"))
    lo,hi=extrema(desired_t2_range); waits=collect(lo:step:hi)
    (waiting_times=waits, phase_cycle=[0,pi/2,pi,3pi/2], coherence_delays=(-0.5,0.5))
end

"""Measure FFT power near `predicted_beat_freq` in a two-column CSV trace."""
function analyze_beats(path::AbstractString; predicted_beat_freq::Real=1.2)
    data=readdlm(path,',',Float64); t=data[:,1]; y=data[:,2].-mean(data[:,2])
    length(t) >= 3 || throw(ArgumentError("at least three samples are required"))
    all(diff(t) .> 0) || throw(ArgumentError("sample times must be strictly increasing"))
    f=FFTW.rfftfreq(length(y),inv(mean(diff(t)))); power=abs2.(FFTW.rfft(y))
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
function _linear_fit_rss(t,y,frequency)
    X = frequency == 0 ? hcat(ones(length(t)),t) : hcat(ones(length(t)),t,sin.(2pi*frequency.*t),cos.(2pi*frequency.*t))
    sum(abs2, y-X*(X\y))
end

"""Analyze a uniformly sampled trace using neutral spectral terminology. Supports
linear detrending and Hann windowing; reports resolution, Nyquist checks, and a
BIC comparison against a no-oscillation linear model. Times and frequencies must
use reciprocal units. A peak is not evidence of quantum coherence by itself."""
function analyze_beats(path::AbstractString; predicted_beat_freq::Real=1.2,
                       detrend::Bool=true, window::Symbol=:hann)
    isfinite(predicted_beat_freq) && predicted_beat_freq >= 0 || throw(ArgumentError("predicted frequency must be finite and nonnegative"))
    data=readdlm(path, ',', Float64); size(data,2)>=2 || throw(ArgumentError("trace must contain two columns"))
    t=data[:,1]; raw=data[:,2]; length(t)>=4 || throw(ArgumentError("at least four samples are required"))
    dt=sampling_interval(t; name="sample times"); all(isfinite,raw) || throw(ArgumentError("signal must be finite"))
    nyquist=1/(2dt); predicted_beat_freq <= nyquist || throw(ArgumentError("predicted frequency exceeds Nyquist frequency $nyquist"))
    duration=t[end]-t[1]; resolution=1/duration
    y=detrend ? raw .- hcat(ones(length(t)),t)*(hcat(ones(length(t)),t)\raw) : raw.-mean(raw)
    w=window == :hann ? 0.5 .- 0.5cos.(2pi*(0:length(y)-1)/(length(y)-1)) : window == :none ? ones(length(y)) : throw(ArgumentError("window must be :hann or :none"))
    f=FFTW.rfftfreq(length(y),inv(dt)); power=abs2.(FFTW.rfft(y.*w)); i=argmin(abs.(f.-predicted_beat_freq))
    rss0=max(_linear_fit_rss(t,raw,0),eps()); rss1=max(_linear_fit_rss(t,raw,f[i]),eps()); n=length(t)
    bic0=n*log(rss0/n)+2log(n); bic1=n*log(rss1/n)+4log(n); delta=bic0-bic1
    warning=duration < 2/max(predicted_beat_freq,eps()) ? "Trace is too short to resolve two predicted beat periods." : nothing
    (; candidate_frequency=f[i], spectral_evidence=(delta_bic=delta, favors_oscillation=delta>0),
       frequency_resolution=resolution, nyquist_frequency=nyquist, duration, warning,
       interpretation="A Fourier candidate is not, by itself, evidence of quantum coherence.")
end

const CODONS=Dict('A'=>["GCT","GCC","GCA","GCG"], 'G'=>["GGT","GGC","GGA","GGG"],
 'L'=>["TTA","TTG","CTT","CTC","CTA","CTG"], 'P'=>["CCT","CCC","CCA","CCG"],
 'F'=>["TTT","TTC"], 'M'=>["ATG"], 'K'=>["AAA","AAG"], 'E'=>["GAA","GAG"],
 'D'=>["GAT","GAC"], 'V'=>["GTT","GTC","GTA","GTG"], 'I'=>["ATT","ATC","ATA"])
const _ADDITIONAL_CODONS = Dict(
    'R'=>["CGT","CGC","CGA","CGG","AGA","AGG"], 'N'=>["AAT","AAC"],
    'C'=>["TGT","TGC"], 'Q'=>["CAA","CAG"], 'H'=>["CAT","CAC"],
    'S'=>["TCT","TCC","TCA","TCG","AGT","AGC"], 'T'=>["ACT","ACC","ACA","ACG"],
    'W'=>["TGG"], 'Y'=>["TAT","TAC"])
merge!(CODONS, _ADDITIONAL_CODONS)
"""Select synonymous codons from a signed relative translation-rate profile."""
function codon_optimizer(sequence::AbstractString, target_profile::AbstractVector{<:Real})
    length(sequence)==length(target_profile) || throw(DimensionMismatch("profile must match sequence"))
    [begin
         choices=get(CODONS,uppercase(aa),nothing)
         isnothing(choices) && throw(ArgumentError("unsupported amino-acid code: $aa"))
 'D'=>["GAT","GAC"], 'V'=>["GTT","GTC","GTA","GTG"], 'I'=>["ATT","ATC","ATA"],
 'R'=>["CGT","CGC","CGA","CGG","AGA","AGG"], 'H'=>["CAT","CAC"],
 'N'=>["AAT","AAC"], 'Q'=>["CAA","CAG"], 'S'=>["TCT","TCC","TCA","TCG","AGT","AGC"],
 'T'=>["ACT","ACC","ACA","ACG"], 'C'=>["TGT","TGC"], 'W'=>["TGG"], 'Y'=>["TAT","TAC"])

"""Select synonymous codons using the sign of a desired local translation profile."""
function codon_optimizer(sequence::AbstractString,target_profile::AbstractVector{<:Real})
    length(sequence)==length(target_profile) || throw(DimensionMismatch("profile must match sequence"))
    [begin
         choices=get(CODONS,uppercase(aa), nothing)
         isnothing(choices) && throw(ArgumentError("unsupported amino acid: $aa"))
         choices[target_profile[i]>0 ? 1 : length(choices)]
     end for (i,aa) in pairs(sequence)]
end
