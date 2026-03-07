# Geant4 Energy Deposition in a Cubic Phantom

This repository contains a **Geant4 simulation** used to obtain 3D energy deposition distributions from **monoenergetic proton beams** in a water phantom.
The simulations serve as the **benchmark reference for the SDE model** used in the accompanying study.

This simulation is part of the MaThRad collaboration work on stochastic differential equation (SDE) models for proton therapy radiation transport.

---

## Overview

- The setup consists of a **20 × 20 × 20 cm³ water phantom** irradiated by a **monoenergetic proton beam**.
- **Energy deposition** is scored as a 3D histogram using a **Custom User Classes** and **G4AnalysisManager**.
- The 3D histogram is created in `RunAction.cc` with a voxel size of 1 cubic mm across the entire phantom.
- The geometry can be set with one of 3 different phantoms: **water, slab** or **insert**, these can be set using UI commands in the macro file.

---

## Physics Lists

To change the physics list, edit `ProtonTherapy.cc` by uncommenting the desired list. The default is QGSP_BIC_EMZ:

```cpp
// Physics lists — uncomment the one you need
auto* physicsFactory = new G4PhysListFactory();

// Reference
auto physicsList = physicsFactory->GetReferencePhysList("QGSP_BIC_EMZ");

// Other physics lists (comment line above and uncomment one of these)
// auto physicsList = physicsFactory->GetReferencePhysList("QGSP_BIC_EMY");
// auto physicsList = physicsFactory->GetReferencePhysList("QGSP_BERT");
```
## Simulation Macros

Each macro modifies the phantom geometry using the command `/phantom/select <geometry>` at the start.
Two proton energies, 100 and 150 MeV, are used in the study. The possible macro configurations are as follows:

| Macro file              | Description                                                       |
|-------------------------|-------------------------------------------------------------------|
| `run0.mac`              | Test run                                                          |
| `run_water_100MeV.mac`  | Pure water phantom, 100 MeV protons                               |
| `run_water_150MeV.mac`  | Pure water phantom, 150 MeV protons                               |
| `run_slab_100MeV.mac`   | 2 cm water + 1 cm bone + 2 cm lung + water, 100 MeV protons       |
| `run_slab_150MeV.mac`   | 2 cm water + 1 cm bone + 2 cm lung + water, 150 MeV protons       |
| `run_insert_100MeV.mac` | 3 cm water + 2 cm bone at top half only + water, 100 MeV protons  |
| `run_insert_150MeV.mac` | 3 cm water + 2 cm bone at top half only + water, 100 MeV protons  |

---

## Compilation

Ensure your Geant4 environment is sourced correctly before compiling:

```bash
source /path/to/geant4-install/bin/geant4.sh
```

Then create a `build/` directory and compile from that folder:
```bash
mkdir build
cd build
cmake ../
make -j$(nproc)
```
## Running the Simulation

After compilation, run the executable by specifying the macro file and number of threads:

```bash
./ProtonTherapy -m <macrofile.mac> -t <nthreads>
```

Each run will generate an output `.root` file containing the **3D energy deposition histogram** in the phantom for post-processing using Python scripts. 
The output files are written to the `Output/` directory by default, but this can be changed in the macro file.

For visualisation only, run the executable without specifying a macro file:
```bash
./ProtonTherapy
```
Then, set the desired geometry with the command `/phantom/select <geometry>` in the Geant4 idle and run the command `/control/execute refresh_vis.mac` to visualise the updated phantom.

## Notes

- This code was developed and tested using **Geant4 version 11.3.0**. Detailed instructions on how to install Geant4 and its prerequisites can be found in https://geant4.web.cern.ch/download/11.3.0.html.
- To record simulation time, run by adding `time` right before the main running command. Example: `time ./ProtonTherapy -m run_water_150MeV.mac -t 8`.