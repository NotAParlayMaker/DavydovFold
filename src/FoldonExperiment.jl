"""Estimate heterodyne SNR with concentration and pulse-energy scaling."""
function predict_signal_amplitude(rnc_concentration::Real, pulse_energy::Real;
                                  efficiency::Real=0.15, noise_floor::Real=1e-6)
    rnc_concentration >= 0 && pulse_energy >= 0 || throw(ArgumentError("inputs must be nonnegative"))
    efficiency >= 0 && noise_floor > 0 || throw(ArgumentError("efficiency must be nonnegative and noise_floor positive"))
    signal=efficiency*rnc_concentration*pulse_energy^1.5
    (signal=signal, snr=signal/noise_floor)
end

"""Design uniformly spaced waiting times and a four-step phase cycle."""
function design_pulse_sequence(desired_t2_range; step::Real=0.05)
    isfinite(step) && step > 0 || throw(ArgumentError("step must be finite and positive"))
    lo,hi=extrema(desired_t2_range); waits=collect(lo:step:hi)
    (waiting_times=waits, phase_cycle=[0,pi/2,pi,3pi/2], coherence_delays=(-0.5,0.5))
end

"""Score a predicted beat in a two-column, uniformly sampled CSV trace."""
function analyze_beats(path::AbstractString; predicted_beat_freq::Real=1.2)
    data=readdlm(path, ',', Float64)
    size(data, 2) >= 2 || throw(ArgumentError("trace must contain at least two columns"))
    t=data[:,1]; y=data[:,2].-mean(data[:,2])
    length(t) >= 4 || throw(ArgumentError("at least four samples are required"))
    all(>(0), diff(t)) || throw(ArgumentError("sample times must be strictly increasing"))
    f=FFTW.rfftfreq(length(y),inv(mean(diff(t)))); power=abs2.(FFTW.rfft(y))
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
