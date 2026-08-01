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
    lo,hi=extrema(desired_t2_range); waits=collect(lo:step:hi)
    (waiting_times=waits, phase_cycle=[0,pi/2,pi,3pi/2], coherence_delays=(-0.5,0.5))
end

"""Measure FFT power near `predicted_beat_freq` in a two-column CSV trace."""
function analyze_beats(path::AbstractString; predicted_beat_freq::Real=1.2)
    data=readdlm(path,',',Float64); t=data[:,1]; y=data[:,2].-mean(data[:,2])
    length(t) >= 3 || throw(ArgumentError("at least three samples are required"))
    all(diff(t) .> 0) || throw(ArgumentError("sample times must be strictly increasing"))
    f=FFTW.rfftfreq(length(y),inv(mean(diff(t)))); power=abs2.(FFTW.rfft(y))
    i=argmin(abs.(f.-predicted_beat_freq)); baseline=median(power[2:end])
    (frequency=f[i], likelihood_ratio=power[i]/max(baseline,eps()), power=power[i])
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
         choices[target_profile[i]>0 ? 1 : length(choices)]
     end for (i,aa) in pairs(sequence)]
end
