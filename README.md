# Proton-beam-SDE
This repository houses a simulator for an SDE describing a the path of a proton in proton beam therapy.
It is under active development and **not fit for standalone use**.

# Installation

In Unix-like environments, simply call `make` in the project root. You will need:

- the g++ compiler (tested on version 13.2.0),
- the [GNU Scientific Library](https://www.gnu.org/software/gsl/) (tested on version 2.7.1),
- the [libconfig](https://hyperrealm.github.io/libconfig/) package (tested on version 1.5).

# Usage

Simulation and output parameters are specified in the config file `sim.cfg`. In particular, you can
modify:
- the number of protons to simulate,
- the radius of the nozzle from which protons are emitted,
- the Gaussian distribution of the initial position of the protons, conditioned to lie in the nozzle,
- the Gaussian distribution of the initial energy of the protons,
- the change points between materials,
- the index of the material in each gap between change points, based on the ordering in Materials/materials.txt,
- the simulator time-step,
- the energy at which a proton is absorbed and its path ends,
- the side length of output voxelation,
- the path to which a sparse 3d grid of dose deposition is stored.

After compilation, run the simulation by calling `./simulate sim.cfg` at the project root.
You can also replace the `sim.cfg` argument with the path to any other appropriate config file.
**Warning**: The default number of protons in `sim.cfg` is set to one million, and voxel sizes and time
steps are also set to a moderate resolution. As a result, the simulation can take 1-2 minutes.

# Post-processing

## Python scripts

Default output paths are specified into the `Output` folder, which contains the Python scripts for output
visualisation. To use it you will need:

- [Python 3](https://www.python.org/) (tested on version 3.12.3), as well as the following Python modules:
- numpy (tested on version 1.26.4),
- matplotlib (tested on version 3.10.0),
- itk (tested on version 5.4.4),
- pymedphys (tested on version 0.40.0),
- scipy (tested on version 1.14.1)
- uproot (tested on version 4.0.0)

All the necessary functions are written in `SDE_vs_G4.py`, whereas the execution scripts to produce the relevant figures 
is written in `SDE_vs_G4_executionScript.py` (if the Geant4 files are available) or `SDE_plots.py` (for exclusively
plotting SDE output). The scripts include integrated and central-axis depth-dose 1D curves, 2D dose distributions, 
lateral profiles, gamma analysis plots, pass rate calculations and proton range calculations. Simply modify the variables 
defined at the start of the script according to what you need to plot. More detailed information is found in the comments
of each script.

## Geant4 output generation

Geant4 is used as reference to assess the accuracy of the SDE. The source code to generate the outputs that are used for 
comparison are found in the `Geant4-dose-and-secondaries` folder, which has its own README file with a detailed explanation
of the model and how to run the simulations.

# Nuclear cross sections

The `Splines` folder contains splines fitted to nuclear cross section data for argon, calcium, carbon, fluorine, hydrogen, nitrogen and oxygen.
When available, the `ENDF/B-VIII.1` nuclear library was used. Otherwise the `JEFF-4.0` library was used. This is sufficient for running simulations
in mediums consisting of air, bone and water. The `Splines` folder also contains three python scripts that are able to generate additional cross sections for input into the code.

- Elastic scattering cross sections for hydrogen can be generated using the script `Elastic_cross_section_generator_hydrogen.py`. The input should be a text file called `Hydrogen_elastic_data.txt` contained within the Splines folder which is in the ENDF format with LAW=5 and LTP=1. The output file is named `hydrogen_el_ruth_cross_sec.txt` which must not be changed so it can be read into the main proton beam code. (Due to heavy recursive calculations this script may take on the order of 10 hours to run.)
- Elastic scattering cross sections for other elements can be generated using the script `Elastic_cross_section_generator.py` The input should be two text files called `Element_cs_el.txt` and `Element_ang_el.txt`. Both files must be contained within the Splines folder and Element is the name of the desired element (the variable Element within the python script must also be replaced with the element's name and all lower case must be used). Both files must be in the ENDF format with LAW=5 and LTP=1, where `Element_cs_el.txt` and `Element_ang_el.txt` correspond to File 3 and File 6 in this format respectively. Furthermore, within the script the variables Z_targ and A_targ must be set equal to the target's charge and mass number respectively. The output file is named `Element_el_ruth_cross_sec.txt` which must not be changed so it can be read into the main proton beam code.
- Inelastic proton-nucleus scattering cross sections can be generated using the script `Non_elastic_cross_section_generator.py`. The input should be two text files called `Element_cs_ne.txt` and `Element_enang_ne.txt`. Both files must be contained within the Splines folder and Element is the name of the desired element (the variable Element within the python script must also be replaced with the element's name and all lower case must be used). Both files must be in the ENDF format using the Kalbach-Mann systematics representation, where `Element_cs_ne.txt` and `Element_enang_ne.txt` correspond to File 3 and File 6 in this format respectively. The output files are named `Element_ne_rate.txt` and `Element_ne_energyangle_cdf.txt` which must not be changed so they can be read into the main proton beam code.

In the case of elastic scattering, the cross section is decomposed into several components which may not all be non-negative. Thus, due to error in the nuclear data, it is possible for the outgoing scattering cdf to be non-monotonically increasing. The python scripts resolve this issue by taking the running supremum of the proposed non-monotonic cdf. When this is required, the code will flag an error to inform the user that this has occured, along with the difference between the non-monotonic cdf and its running supremum.

For those unfamiliar with nuclear data libraries, we give a small tutorial which will allow a user to add additional elements to the proton beam. First navigate to the website `https://www-nds.iaea.org/exfor/endf.htm` where the required nuclear data is housed. 
- For elastic scattering cross sections one should fill in Target: Chemical symbol-Mass number (e.g. O-16 for oxygen 16) Reaction: P,EL, Quantity: SIG for the `Element_cs_el.txt` file and DA/DE for the `Element_ang_el.txt` file. If available, the associated file from the `ENDF/B-VIII.1` nuclear library should be downloaded, otherwise use the `JEFF-4.0` instead.
- For inelastic proton-nucleus scattering cross sections, one should take identical steps except for replacing Reaction: P,EL with Reaction: P,X. Moreover, Quantity: SIG should be used to obtain the `Element_cs_ne.txt` file and DA/DE for `Element_enang_ne.txt` file.

# Materials

The interface for specifying materials is via the text files in the `Materials` folder.

- The `atoms.txt` file lists all chemical elements. Each row specifies an element name (for ease of reference) as well
as its atomic number and mass.
- The `materials.txt` lists names of all desired materials. Each material name must be accompanied by a corresponding
`name.txt` file.
- Each `name.txt` file lists three pieces of information. The first two rows contain the material density in `g / cm^3`
and mean excitation energy in `eV`. Rows from the third specify the chemical composition of the material. Each of these
rows consists of an atomic number and the corresponding fraction by mass of that atom in the material.

For example, water is specified as

1.0\
75.0\
1 0.111\
8 0.889

because it has density 1, mean excitation energy 75, and is approximately 1/9 hydrogen (atomic number 1) and 8/9 oxygen
(atomic number 8) by mass.
