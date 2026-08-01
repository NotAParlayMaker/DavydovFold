#!/usr/bin/env python3
"""Validated, optional OpenMM runner for DavydovFoldon candidate structures."""
from __future__ import annotations
import argparse, json, math, pathlib, sys

REQUIRED = {
    "input_structure_path", "forcefield", "water_model", "temperature_K",
    "pressure_bar", "ionic_strength_molar", "solvent_padding_nm",
    "integration_timestep_fs", "minimization_max_iterations",
    "equilibration_duration_ps", "production_duration_ps",
    "trajectory_output_path", "final_coordinates_path", "run_summary_path", "random_seed",
}

def validate_config(config: dict, config_path: pathlib.Path) -> dict:
    """Validate an OpenMM JSON configuration without importing OpenMM."""
    missing=sorted(REQUIRED-config.keys())
    if missing: raise ValueError("missing required fields: " + ", ".join(missing))
    source=pathlib.Path(config["input_structure_path"])
    if not source.is_absolute(): source=(config_path.parent/source).resolve()
    if not source.is_file(): raise ValueError(f"input structure does not exist: {source}")
    if source.suffix.lower() not in {".pdb", ".cif", ".mmcif"}: raise ValueError("input must be PDB or mmCIF")
    rules={"temperature_K":True,"pressure_bar":False,"ionic_strength_molar":False,
           "solvent_padding_nm":True,"integration_timestep_fs":True,
           "equilibration_duration_ps":False,"production_duration_ps":False}
    for key,positive in rules.items():
        value=config[key]
        if not isinstance(value,(int,float)) or not math.isfinite(value) or (value<=0 if positive else value<0):
            raise ValueError(f"invalid {key}")
    if not isinstance(config["minimization_max_iterations"],int) or config["minimization_max_iterations"]<=0:
        raise ValueError("minimization_max_iterations must be a positive integer")
    config=dict(config); config["input_structure_path"]=str(source); return config

def run(config: dict) -> None:
    """Run minimization, NPT equilibration, and production with optional OpenMM."""
    try:
        from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, unit
        from openmm.app import (DCDReporter, ForceField, HBonds, Modeller, PDBFile,
                                PDBxFile, Simulation, StateDataReporter)
    except ImportError:
        print("OpenMM is not installed. Install it explicitly, e.g. `conda install -c conda-forge openmm`, then rerun.",file=sys.stderr)
        raise SystemExit(3)
    source=pathlib.Path(config["input_structure_path"])
    parser=PDBFile if source.suffix.lower()==".pdb" else PDBxFile
    loaded=parser(str(source)); ff=ForceField(config["forcefield"],config["water_model"])
    modeller=Modeller(loaded.topology,loaded.positions); modeller.addHydrogens(ff)
    modeller.addSolvent(ff,padding=config["solvent_padding_nm"]*unit.nanometer,
        ionicStrength=config["ionic_strength_molar"]*unit.molar)
    system=ff.createSystem(modeller.topology,nonbondedMethod=__import__("openmm.app",fromlist=["PME"]).PME,
        nonbondedCutoff=1.0*unit.nanometer,constraints=HBonds)
    system.addForce(MonteCarloBarostat(config["pressure_bar"]*unit.bar,config["temperature_K"]*unit.kelvin))
    integrator=LangevinMiddleIntegrator(config["temperature_K"]*unit.kelvin,1/unit.picosecond,
        config["integration_timestep_fs"]*unit.femtosecond); integrator.setRandomNumberSeed(config["random_seed"])
    simulation=Simulation(modeller.topology,system,integrator); simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=config["minimization_max_iterations"])
    simulation.context.setVelocitiesToTemperature(config["temperature_K"]*unit.kelvin,config["random_seed"])
    step_fs=config["integration_timestep_fs"]
    equil_steps=round(config["equilibration_duration_ps"]*1000/step_fs)
    prod_steps=round(config["production_duration_ps"]*1000/step_fs); simulation.step(equil_steps)
    interval=max(1,prod_steps//100); simulation.reporters.append(DCDReporter(config["trajectory_output_path"],interval))
    simulation.reporters.append(StateDataReporter(sys.stdout,interval,step=True,potentialEnergy=True,temperature=True))
    simulation.step(prod_steps); state=simulation.context.getState(getPositions=True,getEnergy=True)
    with open(config["final_coordinates_path"],"w") as handle: PDBFile.writeFile(modeller.topology,state.getPositions(),handle)
    summary={"status":"completed","equilibration_steps":equil_steps,"production_steps":prod_steps,
      "potential_energy_kJ_mol":state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
      "scientific_limit":"Short MD and conformation survival do not establish a Davydov mechanism or thermodynamic stability."}
    pathlib.Path(config["run_summary_path"]).write_text(json.dumps(summary,indent=2)+"\n")

def main(argv=None) -> int:
    parser=argparse.ArgumentParser(); parser.add_argument("--config",required=True); parser.add_argument("--validate-only",action="store_true")
    args=parser.parse_args(argv); path=pathlib.Path(args.config)
    try: config=validate_config(json.loads(path.read_text()),path)
    except (OSError,json.JSONDecodeError,ValueError) as exc: print(f"Configuration error: {exc}",file=sys.stderr); return 2
    if args.validate_only: print("Configuration is valid; OpenMM was not imported or executed."); return 0
    run(config); return 0
if __name__=="__main__": raise SystemExit(main())
