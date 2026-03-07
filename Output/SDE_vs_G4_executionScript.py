"""
Calculations and plots for the SDE article
Author: Maria L. Perez-Lara
Returns figures comparing SDE vs Geant4 results in a single specified phantom geometry:
    Figures a and b: 1D slices and projections of the Bragg peaks
    Figures c and d: Lateral profiles at different depths for 100 and 150 MeV
    Figures e and g: 2D slices in z for 100 MeV and 150 MeV
    Figures f and h: Central slice of percentage dose difference array for 100 and 150 MeV
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
    compute_voxelDiff,
)
import numpy as np


"""
Suffix must be selected to compare appropriate files. The options are:
- "_water": For homogeneous water phantom
- "_slab": 2 cm water, 1 cm bone, 2 cm lung, water
- "_insert": 3 cm water, 2 cm bone off axis, water
"""

N_primaries = "1E6"
suffix = "_insert"
energy1 = 100  # MeV
energy2 = 150  # MeV
comparePhysics = False # Set to true if data from BERT and EMY are available to plot comparison between phys lists
voxel_volume = 0.001  # cm3
dd = 2
dta = 0.5
th = 1
fig_extension = ".eps"
MeV_g_to_Gy = 1.60218e-10

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
    material_boundaries = [2, 3, 5]
elif suffix == "_insert":
    density_matrix = include_new_material(
        np.ones((200, 200, 200)), tmin=30, tmax=50, rho=1.45, ycut=100
    )  # add bone
    cuts_100 = [15, 40, 69, 76]
    cuts_150 = [15, 40, 150, 157]
    material_boundaries = [3, 5]
else:
    density_matrix = np.ones((200, 200, 200))
    cuts_100 = [30, 50, 75]
    cuts_150 = [50, 100, 150]
    material_boundaries = None

mass_matrix = density_matrix * voxel_volume
# Data import
g4_dose1 = MeV_g_to_Gy * retrieve_g4_output(f"G4_{N_primaries}_{energy1}MeV{suffix}.root") / mass_matrix
sde_dose1 = MeV_g_to_Gy * retrieve_sde_output(f"SDE_{N_primaries}_{energy1}MeV{suffix}.txt") / mass_matrix

g4_dose2 = MeV_g_to_Gy * retrieve_g4_output(f"G4_{N_primaries}_{energy2}MeV{suffix}.root") / mass_matrix
sde_dose2 = MeV_g_to_Gy * retrieve_sde_output(f"SDE_{N_primaries}_{energy2}MeV{suffix}.txt") / mass_matrix

# Comparative plots - 1D
fa, fb = compare_bragg_peaks(sde_dose1, g4_dose1, energy1, mass_matrix, sde_dose2, g4_dose2, energy2,
                             material_boundaries=material_boundaries)

# Lateral profiles
fc = plot_lateral_profiles(
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

fd = plot_lateral_profiles(
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

# Dose and difference - 2D
fe = plot_slice(
    g4_dose1,
    sde_dose1,
    "z",
    minval=-7,
    xmin=0,
    xmax=9.5,
    ymin=-5,
    ymax=5,
    zmin=-5,
    zmax=5,
)

ff = compute_voxelDiff(
    g4_dose1,
    sde_dose1,
    xmin=0,
    xmax=9.5,
    ymin=-2,
    ymax=2,
    zmin=-5,
    zmax=5,
)

fg = plot_slice(
    g4_dose2,
    sde_dose2,
    "z",
    minval=-7,
    xmin=0,
    xmax=18,
    ymin=-10,
    ymax=10,
    zmin=-10,
    zmax=10,
)

fh = compute_voxelDiff(
    g4_dose2,
    sde_dose2,
    xmin=0,
    xmax=18,
    ymin=-4,
    ymax=4,
    zmin=-10,
    zmax=10,
)

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

# Save figures with corresponding names according to PMB article
if suffix == "_water":
    fa.savefig(f"figure2a{fig_extension}")
    fb.savefig(f"figure2b{fig_extension}")
    fc.savefig(f"figure3a{fig_extension}")
    fd.savefig(f"figure3b{fig_extension}")
    fe.savefig(f"figure4a{fig_extension}")
    ff.savefig(f"figure4b{fig_extension}")
    fg.savefig(f"figure4c{fig_extension}")
    fh.savefig(f"figure4d{fig_extension}")
if suffix == "_slab":
    fa.savefig(f"figure6a{fig_extension}")
    fb.savefig(f"figure6b{fig_extension}")
    fc.savefig(f"figure7a{fig_extension}")
    fd.savefig(f"figure7b{fig_extension}")
    fe.savefig(f"figure8a{fig_extension}")
    ff.savefig(f"figure8b{fig_extension}")
    fg.savefig(f"figure8c{fig_extension}")
    fh.savefig(f"figure8d{fig_extension}")
if suffix == "_insert":
    fa.savefig(f"figure9a{fig_extension}")
    fb.savefig(f"figure9b{fig_extension}")
    fc.savefig(f"figure10a{fig_extension}")
    fd.savefig(f"figure10b{fig_extension}")
    fe.savefig(f"figure11a{fig_extension}")
    ff.savefig(f"figure11b{fig_extension}")
    fg.savefig(f"figure11c{fig_extension}")
    fh.savefig(f"figure11d{fig_extension}")

# Comparison with other physics lists, only available for homogeneous phantom configs:
if comparePhysics and suffix == "_water":
    g4_dose_emy1 = MeV_g_to_Gy * retrieve_g4_output(
        f"G4_{N_primaries}_{energy1}MeV{suffix}_emy.root"
    ) / mass_matrix
    g4_dose_bert1 = MeV_g_to_Gy * retrieve_g4_output(
        f"G4_{N_primaries}_{energy1}MeV{suffix}_bert.root"
    ) / mass_matrix
    g4_dose_emy2 = MeV_g_to_Gy * retrieve_g4_output(
        f"G4_{N_primaries}_{energy2}MeV{suffix}_emy.root"
    ) / mass_matrix
    g4_dose_bert2 = MeV_g_to_Gy * retrieve_g4_output(
        f"G4_{N_primaries}_{energy2}MeV{suffix}_bert.root"
    ) / mass_matrix

    fi = plot_multiple_bragg_peaks(
        g4_dose1,
        g4_dose_emy1,
        g4_dose_bert1,
        sde_dose1,
        mass_matrix=mass_matrix,
        how="proj",
        xmax=8,
        names=["QGSP_BIC_EMY", "QGSP_BERT", "SDE"],
        ref_name="QGSP_BIC_EMZ",
        maxdif=10,
    )
    fj = plot_multiple_bragg_peaks(
        g4_dose2,
        g4_dose_emy2,
        g4_dose_bert2,
        sde_dose2,
        mass_matrix=mass_matrix,
        how="proj",
        xmax=17,
        names=["QGSP_BIC_EMY", "QGSP_BERT", "SDE"],
        ref_name="QGSP_BIC_EMZ",
        maxdif=25,
    )

    fi.savefig(f"figure5a{fig_extension}")
    fj.savefig(f"figure5b{fig_extension}")

    # Uncomment the following lines to plot 1D slices
    # fk = plot_multiple_bragg_peaks(
    #     g4_dose1,
    #     g4_dose_emy1,
    #     g4_dose_bert1,
    #     sde_dose1,
    #     mass_matrix=mass_matrix,
    #     how="slice",
    #     xmax=8,
    #     names=["QGSP_BIC_EMY", "QGSP_BERT", "SDE"],
    #     ref_name="QGSP_BIC_EMZ",
    # )
    # fl = plot_multiple_bragg_peaks(
    #     g4_dose2,
    #     g4_dose_emy2,
    #     g4_dose_bert2,
    #     sde_dose2,
    #     mass_matrix=mass_matrix,
    #     how="slice",
    #     xmax=17,
    #     names=["QGSP_BIC_EMY", "QGSP_BERT", "SDE"],
    #     ref_name="QGSP_BIC_EMZ",
    #     maxdif=30,
    # )

