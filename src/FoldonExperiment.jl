"""Estimate heterodyne SNR with concentration and pulse-energy scaling."""
function predict_signal_amplitude(rnc_concentration, pulse_energy; efficiency=0.15, noise_floor=1e-6)
    rnc_concentration >= 0 && pulse_energy >= 0 || throw(ArgumentError("inputs must be nonnegative"))
    signal=efficiency*rnc_concentration*pulse_energy^1.5
    (signal=signal, snr=signal/noise_floor)
end

function design_pulse_sequence(desired_t2_range; step=0.05)
    lo,hi=extrema(desired_t2_range); waits=collect(lo:step:hi)
    (waiting_times=waits, phase_cycle=[0,pi/2,pi,3pi/2], coherence_delays=(-0.5,0.5))
end

function analyze_beats(path::AbstractString; predicted_beat_freq=1.2)
    data=readdlm(path,',',Float64); t=data[:,1]; y=data[:,2].-mean(data[:,2])
    f=FFTW.rfftfreq(length(y),inv(mean(diff(t)))); power=abs2.(FFTW.rfft(y))
    i=argmin(abs.(f.-predicted_beat_freq)); baseline=median(power[2:end])
    (frequency=f[i], likelihood_ratio=power[i]/max(baseline,eps()), power=power[i])
end

const CODONS=Dict('A'=>["GCT","GCC","GCA","GCG"], 'G'=>["GGT","GGC","GGA","GGG"],
 'L'=>["TTA","TTG","CTT","CTC","CTA","CTG"], 'P'=>["CCT","CCC","CCA","CCG"],
 'F'=>["TTT","TTC"], 'M'=>["ATG"], 'K'=>["AAA","AAG"], 'E'=>["GAA","GAG"],
 'D'=>["GAT","GAC"], 'V'=>["GTT","GTC","GTA","GTG"], 'I'=>["ATT","ATC","ATA"])
function codon_optimizer(sequence::AbstractString,target_profile)
    length(sequence)==length(target_profile) || throw(DimensionMismatch("profile must match sequence"))
    [begin
         choices=get(CODONS,uppercase(aa),["NNN"])
         choices[target_profile[i]>0 ? 1 : length(choices)]
     end for (i,aa) in pairs(sequence)]
end
