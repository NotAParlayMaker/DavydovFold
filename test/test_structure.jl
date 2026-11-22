@testset "candidate backbone reconstruction" begin
    one=reconstruct_backbone("G",[-60.0],[-45.0])
    @test length(one.atoms)==4
    @test one.diagnostics.radius_of_gyration_A==0.0
    @test all(a->all(isfinite,a.xyz_A),one.atoms)

    peptide=reconstruct_backbone("AGST",fill(-60.0,4),fill(-45.0,4))
    again=reconstruct_backbone("AGST",fill(-60.0,4),fill(-45.0,4))
    @test length(peptide.atoms)==16
    @test getfield.(peptide.atoms,:xyz_A)==getfield.(again.atoms,:xyz_A)
    @test size(peptide.diagnostics.contact_map)==(4,4)
    @test length(peptide.diagnostics.secondary_structure)==4
    @test_throws DimensionMismatch reconstruct_backbone("AG",[-60.0],[-45.0,-45.0])
    @test_throws ArgumentError reconstruct_backbone("AG",[-Inf,-60.0],[-45.0,-45.0])
    @test_throws ArgumentError reconstruct_backbone("AG",[-361.0,-60.0],[-45.0,-45.0])

    mktempdir() do dir
        pdb=joinpath(dir,"candidate.pdb"); cif=joinpath(dir,"candidate.cif")
        @test write_pdb(pdb,peptide)==pdb; @test write_mmcif(cif,peptide)==cif
        lines=readlines(pdb); @test startswith(lines[3],"ATOM  ")
        @test any(contains("MODEL-DERIVED"),lines); @test lines[end]=="END"
        @test contains(read(cif,String),"_atom_site.Cartn_x")
    end

    # A deliberately folded idealized chain exercises deterministic clash reporting.
    compact=reconstruct_backbone("AAAAAAAA",fill(-180.0,8),fill(0.0,8);clash_cutoff_A=3.0)
    @test !isempty(compact.diagnostics.clashes)
end

@testset "explicit reduced torsion mapping" begin
    traj,_,_=run_soliton("AGS";tmax=0.0)
    mapped=approximate_backbone_torsions(traj,:last)
    @test length(mapped.torsions.phi)==3
    @test !isempty(mapped.warnings)
    structure=reconstruct_structure(traj)
    @test structure.metadata["source"]=="reduced_trajectory_approximate_mapping"
    @test any(contains("does not uniquely determine"),structure.diagnostics.warnings)
end
