@testset "OpenMM export specification" begin
    mktempdir() do dir
        pdb=joinpath(dir,"candidate.pdb")
        structure=reconstruct_backbone("AG",fill(-60.0,2),fill(-45.0,2)); write_pdb(pdb,structure)
        spec=prepare_openmm_system(structure;input_structure_path=pdb,production_duration_ps=1.0)
        out=joinpath(dir,"system.json"); @test write_openmm_spec(out,spec)==out
        text=read(out,String)
        for field in ("temperature_K","pressure_bar","integration_timestep_fs","production_duration_ps")
            @test contains(text,field)
        end
        @test contains(text,"Short MD does not prove thermodynamic stability")
        @test_throws ArgumentError prepare_openmm_system(pdb;temperature_K=0.0)
        @test_throws ArgumentError prepare_openmm_system(joinpath(dir,"missing.pdb"))
    end
end
