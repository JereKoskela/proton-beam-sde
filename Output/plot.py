import numpy as np
import matplotlib.pyplot as plt
import itk


# Obtain SDE output from the txt file generated
def retrieve_sde_output(filename, airgap=0.1, bin_size=0.1):
    data = np.genfromtxt(filename)
    # Filter to region of interest, adjusted by airgap on x
    data = data[
        (data[:, 0] >= airgap)
        & (data[:, 0] < airgap + 10)
        & (data[:, 1] >= -5)
        & (data[:, 1] < 5)
        & (data[:, 2] >= -5)
        & (data[:, 2] < 5)
    ]
    x, y, z, vals = data.T
    x = x - airgap
    x_bins = int(10 / bin_size)
    y_bins = int(10 / bin_size)
    z_bins = int(10 / bin_size)
    volume = np.zeros((x_bins, y_bins, z_bins))

    for i in range(len(x)):
        xi = int(np.round(x[i] / bin_size))
        yi = int(np.round((y[i] + 5) / bin_size))  # shift y to [0,10]
        zi = int(np.round((z[i] + 5) / bin_size))  # shift z to [0,10]
        volume[xi, yi, zi] += vals[i]

    return volume


# Compares 2D projections
def plot_projection(sde_mat):
    sde_dose_proj = np.sum(sde_mat, axis=2)

    # For a common colorbar ignoring infinite values when taking the log
    with np.errstate(divide="ignore"):
        sde_log = np.log10(sde_dose_proj)

    sde_log[np.isneginf(sde_log)] = np.nan

    vmin = np.nanmin(sde_log)
    vmax = np.nanmax(sde_log)

    plt.plot()
    im = plt.imshow(
        sde_log.T, extent=[0, 10, -5, 5], origin="lower", vmin=vmin, vmax=vmax
    )
    plt.colorbar(
        im, label=r"Log(Dose[$\mu$Gy]), summed over z", fraction=0.05, pad=0.04
    )
    plt.ylabel("y (cm)")
    plt.xlabel("x (cm)")
    plt.title("X-Y projection", fontweight="bold")
    plt.show(block=True)


# Compares 1D dose either summed (projection) or from a slice
def plot_bragg_peak(
    sde_mat, how="proj"
):  # Define if projected Bragg peak or from a slice
    if how == "proj":
        sde_dose1d = np.sum(sde_mat, axis=(1, 2))
    else:
        sde_dose1d = sde_mat[:, 49, 49]

    x = np.linspace(0.1, 10, 100)
    plt.figure()
    plt.plot(x, sde_dose1d, label="SDE")
    plt.xlabel("Depth (cm)")
    plt.ylabel(r"Dose ($\mu$Gy)")
    plt.show(block=True)


# Shows all central slices for one 3D dose distribution only
def plot_all_slices(mat):
    slicex = mat[49, :, :]
    slicey = mat[:, 49, :]
    slicez = mat[:, :, 49]
    # For a common colorbar ignoring infinite values when taking the log
    with np.errstate(divide="ignore"):
        slicex_log = np.log10(slicex)
        slicey_log = np.log10(slicey)
        slicez_log = np.log10(slicez)

    slicex_log[np.isneginf(slicex_log)] = np.nan
    slicey_log[np.isneginf(slicey_log)] = np.nan
    slicez_log[np.isneginf(slicez_log)] = np.nan

    vmin = np.nanmin([slicex_log, slicey_log, slicez_log])
    vmax = np.nanmax([slicex_log, slicey_log, slicez_log])

    f, ax = plt.subplots(1, 3, sharey=True, figsize=(20, 7))
    im1 = ax[0].imshow(
        slicex_log.T, extent=[-5, 5, -5, 5], origin="lower", vmin=vmin, vmax=vmax
    )
    im2 = ax[1].imshow(
        slicey_log.T, extent=[0, 10, -5, 5], origin="lower", vmin=vmin, vmax=vmax
    )
    im3 = ax[2].imshow(
        slicez_log.T, extent=[0, 10, -5, 5], origin="lower", vmin=vmin, vmax=vmax
    )
    f.colorbar(
        im3,
        ax=ax.ravel().tolist(),
        label=r"Log(Dose[$\mu$Gy])",
        fraction=0.01,
        pad=0.04,
    )
    ax[0].set_title("X central slice", fontweight="bold")
    ax[1].set_title("Y central slice", fontweight="bold")
    ax[2].set_title("Z central slice", fontweight="bold")
    plt.subplots_adjust(
        top=1, bottom=0.1, right=0.88, left=0.05, hspace=0.05, wspace=0.1
    )
    plt.show(block=True)


# ---------------------------TODO: EXECUTION --------------------------------------
MeV_to_Gy = 1.602 * 10 ** (-13)

sde = retrieve_sde_output("Output/test.txt")
sde_dose = sde * MeV_to_Gy * 10**6  # conversion to dose in microGy

plot_projection(sde_dose)
plot_all_slices(sde_dose)
plot_bragg_peak(sde_dose, how="proj")
plot_bragg_peak(sde_dose, how="slice")
