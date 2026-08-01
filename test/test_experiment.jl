@testset "experiment design" begin
    @test predict_signal_amplitude(2.0, 3.0).snr > 0
    @test length(design_pulse_sequence((0.0, 1.0); step=0.1).waiting_times) == 11
    @test codon_optimizer("MAG", [1, -1, 1]) == ["ATG", "GCG", "GGT"]
    @test_throws DimensionMismatch codon_optimizer("MA", [1])
    @test codon_optimizer("WR", [1, -1]) == ["TGG", "AGG"]
    @test_throws ArgumentError codon_optimizer("X", [1])
    @test_throws ArgumentError design_pulse_sequence((0.0, 1.0); step=Inf)

    mktemp() do path, io
        writedlm(io, [0.0 0.0; 0.1 1.0; 0.25 0.0; 0.3 -1.0], ',')
        close(io)
        @test_throws ArgumentError analyze_beats(path)
    end
    @test_throws ArgumentError analyze_beats("unused.csv"; predicted_beat_freq=Inf)
end
