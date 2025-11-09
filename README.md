# Proton-beam-SDE
This repository houses a simulator for an SDE describing a the path of a proton in proton beam therapy.
This branch is records the version used for the simulations in
[Crossley et al. (2025)](https://doi.org/10.1098/rspa.2024.0687). See the main branch for an
up-to-date implementation.

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
- the distribution of the initial energy of the protons,
- the simulator time-step,
- the minimum elastic scatter angle (scattering at lower angles being subsumed into a diffusion term),
- the energy at which a proton is absorbed and its path ends,
- the side length of output voxelation (**Warning**: very small voxels result in a large memory cost,
  especially if 3d output is desired),
- desired output modes:
  * `output_3d = 1` outputs the full 3d voxelation of dose deposition,
  * `output_2d = 1` and `output_1d = 1` output the dose with 1 or 2 dimensions integrated out,
  * `output_2d_slice = 1` returns the 2d dose deposition along a one-voxel slice at `z = 0`.
  * Setting any output mode to 0 means that output is not produced. A 3d voxel grid is only created if
    `output_3d = 1`, which slows down the simulation and increases memory cost significantly.
  * The `out_path` variables specify files in which the outputs should be stored.

After compilation, run the simulation by calling `./simulate sim.cfg` at the project root.
You can also replace the `sim.cfg` argument with the path to any other appropriate config file.
**Warning**: The default number of protons in `sim.cfg` is set to one million, and voxel sizes and time
steps are also set to a rather high resolution. As a result, the simulation can take around an hour to run.

# Post-processing

Visualisation tools for `output_2d` `output_1d`, and `output_2d_slice` are provided in `plots.R`. To use
them you will need:

- [R](https://www.r-project.org/) (tested on version 4.3.3),
- the [viridis](https://cran.r-project.org/web/packages/viridis/index.html) package (tested on version 0.6.5),
- while not strictly necessary, [RStudio](https://posit.co/) or a similar IDE will also make life easier
  (tested on RStudio version 2023.06.1).

If you have used the default number of protons, as well as output modes and paths in `sim.cfg`, you can simply
run `plots.R` to produce heatmaps of the integrated 2d dose, the 2d dose slice, and the integrated 1d dose.
Make sure the R working directory is set to the project root (or the location of your simulation output files,
if different).
