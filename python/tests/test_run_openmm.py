import importlib.util, json, pathlib, subprocess, sys

ROOT=pathlib.Path(__file__).parents[2]
SPEC=importlib.util.spec_from_file_location("runner",ROOT/"python/run_openmm.py")
runner=importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(runner)

def config(path):
    return {"input_structure_path":str(path),"forcefield":"amber14-all.xml","water_model":"amber14/tip3p.xml",
      "temperature_K":310.0,"pressure_bar":1.0,"ionic_strength_molar":0.15,"solvent_padding_nm":1.0,
      "integration_timestep_fs":2.0,"minimization_max_iterations":100,"equilibration_duration_ps":1.0,
      "production_duration_ps":1.0,"trajectory_output_path":"trajectory.dcd","final_coordinates_path":"final.pdb",
      "run_summary_path":"summary.json","random_seed":42}

def test_validation_does_not_require_openmm(tmp_path):
    pdb=tmp_path/"candidate.pdb"; pdb.write_text("END\n")
    cfg=tmp_path/"config.json"; cfg.write_text(json.dumps(config(pdb)))
    result=subprocess.run([sys.executable,str(ROOT/"python/run_openmm.py"),"--config",str(cfg),"--validate-only"],capture_output=True,text=True)
    assert result.returncode==0 and "not imported" in result.stdout

def test_rejects_missing_structure(tmp_path):
    cfg=config(tmp_path/"missing.pdb")
    try: runner.validate_config(cfg,tmp_path/"config.json")
    except ValueError as exc: assert "does not exist" in str(exc)
    else: raise AssertionError("missing structure accepted")
