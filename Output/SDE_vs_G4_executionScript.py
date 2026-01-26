"""
Calculations and plots for the SDE article
Author: Maria L. Perez-Lara
Returns figures comparing SDE vs Geant4 results in a single specified phantom geometry:
    Figures 1 and 4: 2D slices in z for 100 MeV and 150 MeV
    Figures 2 and 5: Lateral profiles at different depths for 100 and 150 MeV
    Figures 3 and 6: Central slice of percentage dose difference array for 100 and 150 MeV
    Figures 7 and 8: 1D slices and projections of the Bragg peaks
    Other figures: In case comparePhysics is set to True
Also prints in the terminal the following values:
    Gamma pass rates using the DD, DTA and TH values specified
    Range agreement for 100 MeV and 150 MeV
"""

from SDE_vs_G4 import (
    retrieve_sde_output,
    retrieve_g4_output,
    plot_lateral_profiles,
    plot_slice,
    pymedphys_gamma,
    compare_bragg_peaks,
    include_new_material,
    range_comparison,
    plot_multiple_bragg_peaks,
    compute_percentage_voxelDiff,
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
energy1 = 100  # MeV
energy2 = 150  # MeV
comparePhysics = False  # Set to true if data from BERT and EMY are available to plot comparison between phys lists
voxel_volume = 0.001  # cm3
dd = 2
dta = 0.5
th = 1

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

# Data import
g4_dose1 = retrieve_g4_output(f"G4_{N_primaries}_{energy1}MeV{suffix}.root") / (
    float(N_primaries) * density_matrix * voxel_volume
)
sde_dose1 = retrieve_sde_output(f"SDE_{N_primaries}_{energy1}MeV{suffix}.txt") / (
    float(N_primaries) * density_matrix * voxel_volume
)

g4_dose2 = retrieve_g4_output(f"G4_{N_primaries}_{energy2}MeV{suffix}.root") / (
    float(N_primaries) * density_matrix * voxel_volume
)
sde_dose2 = retrieve_sde_output(f"SDE_{N_primaries}_{energy2}MeV{suffix}.txt") / (
    float(N_primaries) * density_matrix * voxel_volume
)

# Individual energy plots - energy 1
f1 = plot_slice(
    g4_dose1,
    sde_dose1,
    "z",
    minval=-5,
    xmin=0,
    xmax=9.5,
    ymin=-5,
    ymax=5,
    zmin=-5,
    zmax=5,
)
f1.savefig(f"2DSlice_z_{energy1}MeV{suffix}.png")

f2 = plot_lateral_profiles(
    g4_dose1,
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
f2.savefig(f"lateral_profiles_{energy1}MeV{suffix}.png")

f3 = compute_percentage_voxelDiff(
    g4_dose1,
    sde_dose1,
    xmin=0,
    xmax=9.5,
    ymin=-2,
    ymax=2,
    zmin=-5,
    zmax=5,
)
f3.savefig(f"doseDiff_{energy1}MeV{suffix}.png")

# Gamma analysis
print("Gamma values for 100 MeV: ")
g1, pr1 = pymedphys_gamma(
    g4_dose1,
    sde_dose1,
    dta=dta,
    dd=dd,
    th_percent=th,
    is_local=True,
    xmin=0,
    xmax=9.5,
    ymin=-2,
    ymax=2,
    zmin=-5,
    zmax=5,
)

# Individual plots - energy 2
f4 = plot_slice(
    g4_dose2,
    sde_dose2,
    "z",
    minval=-5,
    xmin=0,
    xmax=18,
    ymin=-10,
    ymax=10,
    zmin=-10,
    zmax=10,
)
f4.savefig(f"2DSlice_z_{energy2}MeV{suffix}.png")

f5 = plot_lateral_profiles(
    g4_dose2,
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
f5.savefig(f"lateral_profiles_{energy2}MeV{suffix}.png")

f6 = compute_percentage_voxelDiff(
    g4_dose2,
    sde_dose2,
    xmin=0,
    xmax=18,
    ymin=-4,
    ymax=4,
    zmin=-10,
    zmax=10,
)
f6.savefig(f"doseDiff_{energy2}MeV{suffix}.png")

# Gamma analysis
print("Gamma values for 150 MeV: ")
g2, pr2 = pymedphys_gamma(
    g4_dose2,
    sde_dose2,
    dta=dta,
    dd=dd,
    th_percent=th,
    is_local=True,
    xmin=0,
    xmax=18,
    ymin=-4,
    ymax=4,
    zmin=-5,
    zmax=5,
)

# Comparative plots
f7, f8 = compare_bragg_peaks(sde_dose1, g4_dose1, energy1, sde_dose2, g4_dose2, energy2)
f7.savefig(f"1DProj_{energy1}_vs_{energy2}MeV{suffix}.png")
f8.savefig(f"1DSlice_{energy1}_vs_{energy2}MeV{suffix}.png")

# Proton range numbers
r90_1 = range_comparison(g4_dose1, sde_dose1, 0.9)
r50_1 = range_comparison(g4_dose1, sde_dose1, 0.5)
print(
    f"For {energy1} MeV, proton range agreement is within {r90_1*10:.2f} mm (R90) and {r50_1*10:.2f} mm (R50)"
)
r90_2 = range_comparison(g4_dose2, sde_dose2, 0.9)
r50_2 = range_comparison(g4_dose2, sde_dose2, 0.5)
print(
    f"For {energy2} MeV, proton range agreement is within {r90_2*10:.2f} mm (R90) and {r50_2*10:.2f} mm (R50)"
)

# Comparison with other physics lists, only available for homogeneous phantom configs:
if comparePhysics and suffix == "_water":
    g4_dose_emy1 = retrieve_g4_output(
        f"G4_{N_primaries}_{energy1}MeV{suffix}_emy.root"
    ) / (float(N_primaries) * density_matrix * voxel_volume)
    g4_dose_bert1 = retrieve_g4_output(
        f"G4_{N_primaries}_{energy1}MeV{suffix}_bert.root"
    ) / (float(N_primaries) * density_matrix * voxel_volume)
    g4_dose_emy2 = retrieve_g4_output(
        f"G4_{N_primaries}_{energy2}MeV{suffix}_emy.root"
    ) / (float(N_primaries) * density_matrix * voxel_volume)
    g4_dose_bert2 = retrieve_g4_output(
        f"G4_{N_primaries}_{energy2}MeV{suffix}_bert.root"
    ) / (float(N_primaries) * density_matrix * voxel_volume)
    f9 = plot_multiple_bragg_peaks(
        g4_dose1,
        g4_dose_emy1,
        g4_dose_bert1,
        sde_dose1,
        how="proj",
        xmax=8,
        names=["QGSP_BIC_EMY", "QGSP_BERT", "SDE"],
        ref_name="QGSP_BIC_EMZ",
        maxdif=10,
    )
    f9.savefig(f"1DProj_multiComparison_{energy1}MeV.png")

    f10 = plot_multiple_bragg_peaks(
        g4_dose1,
        g4_dose_emy1,
        g4_dose_bert1,
        sde_dose1,
        how="slice",
        xmax=8,
        names=["QGSP_BIC_EMY", "QGSP_BERT", "SDE"],
        ref_name="QGSP_BIC_EMZ",
    )
    f10.savefig(f"1DSlice_multiComparison_{energy1}MeV.png")

    f11 = plot_multiple_bragg_peaks(
        g4_dose2,
        g4_dose_emy2,
        g4_dose_bert2,
        sde_dose2,
        how="proj",
        xmax=17,
        names=["QGSP_BIC_EMY", "QGSP_BERT", "SDE"],
        ref_name="QGSP_BIC_EMZ",
        maxdif=25,
    )
    f11.savefig(f"1DProj_multiComparison_{energy2}MeV.png")

    f12 = plot_multiple_bragg_peaks(
        g4_dose2,
        g4_dose_emy2,
        g4_dose_bert2,
        sde_dose2,
        how="slice",
        xmax=17,
        names=["QGSP_BIC_EMY", "QGSP_BERT", "SDE"],
        ref_name="QGSP_BIC_EMZ",
        maxdif=30,
    )
    f12.savefig(f"1DSlice_multiComparison_{energy2}MeV.png")
