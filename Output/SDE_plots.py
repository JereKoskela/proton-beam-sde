"""
Calculations and plots for the SDE article
This script does not include Geant4 data, it is used for SDE data visualisation purposes.
Returns:
    Figure with 2D slices in x, y, z
    Figure with lateral profiles at different depths
    Figure with 1D projections of the bragg peaks
    Figure with 1D slices of the bragg peaks
Author: Maria L. Perez-Lara
"""

from SDE_vs_G4 import (
    retrieve_sde_output,
    plot_lateral_profiles,
    include_new_material,
    plot_bragg_peaks_SDE,
    plot_all_slices,
)
import numpy as np


"""
Suffix must be selected to compare appropriate files. The options are:
- "_water": For homogeneous water phantom
- "_slab": 2 cm water, 1 cm bone, 2 cm lung, water
- "_insert": 3 cm water, 2 cm bone off axis, water
"""

N_primaries = "1E6"
suffix = "_water"
energy1 = 100  # MeV, must always have a value
energy2 = None  # MeV, can be set to None to skip. Only energy1 plots will be generated.
voxel_volume = 0.001  # cm3

# Adjust parameters according to case:
if suffix == "_slab":
    density_matrix = include_new_material(
        np.ones((200, 200, 200)), tmin=20, tmax=30, rho=1.45
    )  # add bone
    density_matrix = include_new_material(
        density_matrix, tmin=30, tmax=50, rho=0.385
    )  # add lung
    cuts_100 = [10, 25, 40, 85]
    cuts_150 = [10, 25, 40, 165]
elif suffix == "_insert":
    density_matrix = include_new_material(
        np.ones((200, 200, 200)), tmin=30, tmax=50, rho=1.45, ycut=100
    )  # add bone
    cuts_100 = [15, 40, 69, 76]
    cuts_150 = [15, 40, 150, 157]
else:
    density_matrix = np.ones((200, 200, 200))
    cuts_100 = [30, 50, 75]
    cuts_150 = [50, 100, 150]

sde_dose1 = retrieve_sde_output(f"SDE_{N_primaries}_{energy1}MeV{suffix}.txt") / (
    float(N_primaries) * density_matrix * voxel_volume
)
f1 = plot_all_slices(
    sde_dose1,
    xmin=0,
    xmax=9.5,
    ymin=-5,
    ymax=5,
    zmin=-5,
    zmax=5,
)
f1.savefig(f"../Output/SDE_2DSlices_{energy1}MeV{suffix}.png")

f2 = plot_lateral_profiles(
    None,
    sde_dose1,
    xmin=0,
    xmax=10,
    ymin=-5,
    ymax=5,
    zmin=-5,
    zmax=5,
    lowlim=-1,
    uplim=1,
    cuts=cuts_100,
)
f2.savefig(f"../Output/SDE_lateral_profiles_{energy1}MeV{suffix}.png")

if energy2 is None:
    f3, f4 = plot_bragg_peaks_SDE(sde_dose1, energy1)
    f3.savefig(f"../Output/SDE_1DProj_{energy1}MeV{suffix}.png")
    f4.savefig(f"../Output/SDE_1DSlice_{energy1}MeV{suffix}.png")

else:
    sde_dose2 = retrieve_sde_output(f"SDE_{N_primaries}_{energy2}MeV{suffix}.txt") / (
        float(N_primaries) * density_matrix * voxel_volume
    )
    f5 = plot_all_slices(
        sde_dose2,
        xmin=0,
        xmax=18,
        ymin=-10,
        ymax=10,
        zmin=-10,
        zmax=10,
    )
    f5.savefig(f"../Output/SDE_2DSlice_z_{energy2}MeV{suffix}.png")

    f6 = plot_lateral_profiles(
        None,
        sde_dose2,
        xmin=0,
        xmax=18,
        ymin=-10,
        ymax=10,
        zmin=-10,
        zmax=10,
        lowlim=-1,
        uplim=1,
        cuts=cuts_150,
    )
    f6.savefig(f"../Output/SDE_lateral_profiles_{energy2}MeV{suffix}.png")

    f7, f8 = plot_bragg_peaks_SDE(sde_dose1, energy1, sde_dose2, energy2)
    f7.savefig(f"../Output/SDE_1DProj_{energy1}_vs_{energy2}MeV{suffix}.png")
    f8.savefig(f"../Output/SDE_1DSlice_{energy1}_vs{energy2}MeV{suffix}.png")
