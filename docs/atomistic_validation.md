# Atomistic export and validation boundary

DavydovFoldon writes a unit-labeled JSON specification for an optional external
OpenMM process. It neither imports OpenMM during Julia package loading nor
installs or launches it automatically. Use:

```bash
python python/run_openmm.py --config examples/openmm_system.json --validate-only
python python/run_openmm.py --config system.json
```

The runner adds supported hydrogens, explicit solvent and ions, minimizes,
equilibrates, runs production, and writes DCD, final PDB, and JSON summary files.
A force field is an approximation. Short MD does not prove thermodynamic
stability; survival does not establish a Davydov mechanism; failure can reflect
idealized reconstruction or force-field limitations.

Atomistic trajectory import and quantitative comparison remain follow-up work.
