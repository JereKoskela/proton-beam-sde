"""
Functions for the SDE article
Author: Maria L. Perez-Lara
"""

import numpy as np
import matplotlib.pyplot as plt
import pymedphys
from scipy.interpolate import interp1d
from matplotlib.lines import Line2D
import uproot4

plt.rcParams.update({"font.size": 13})


# Obtain SDE output from the txt file generated. Based on 20x20x20 cm3 cubic phantom area.
# Set min and max values for each axis in cm
def retrieve_sde_output(
    filename, bin_size=0.1, xmin=0, xmax=20, ymin=-10, ymax=10, zmin=-10, zmax=10
):
    data = np.genfromtxt(filename)
    # Filter to region of interest (20x20x20 cm3 phantom)
    data = data[
        (data[:, 0] >= xmin)
        & (data[:, 0] < xmax)
        & (data[:, 1] >= ymin)
        & (data[:, 1] < ymax)
        & (data[:, 2] >= zmin)
        & (data[:, 2] < zmax)
    ]
    x, y, z, vals = data.T
    x_bins = int((xmax - xmin) / bin_size)
    y_bins = int((xmax - xmin) / bin_size)
    z_bins = int((xmax - xmin) / bin_size)
    volume = np.zeros((x_bins, y_bins, z_bins))

    for i in range(len(x)):
        xi = int(np.round(x[i] / bin_size))
        yi = int(
            np.round((y[i] + int((ymax - ymin) / 2)) / bin_size)
        )  # shift y to [0,20]
        zi = int(
            np.round((z[i] + int(zmax - zmin) / 2) / bin_size)
        )  # shift z to [0,20]
        volume[xi, yi, zi] += vals[i]
    return volume


# Obtain geant4 output from the root file generated. Based on 20x20x20 cm3 cubic phantom by default (voxel size 1 mm)
def retrieve_g4_output(filename):
    data = uproot4.open(filename)
    hist = data["TotalEneDep"]
    values = hist.values()
    return values


# Cut down a 3D array to the region of interest
def define_ROI(
    mat,
    xmin=0,
    xmax=20,
    ymin=-10,
    ymax=10,
    zmin=-10,
    zmax=10,
    voxelOrigin=np.array([0, -10, -10]),
):
    minArray = np.array([xmin, ymin, zmin])
    maxArray = np.array([xmax, ymax, zmax])
    minlims = (minArray - voxelOrigin) / 0.1
    maxlims = (maxArray - voxelOrigin) / 0.1
    mat = mat[
        int(minlims[0]) : int(maxlims[0]),
        int(minlims[1]) : int(maxlims[1]),
        int(minlims[2]) : int(maxlims[2]),
    ]
    return mat


# Compares central slice in a defined axis for both geant4 and SDE
# set a minval to change the minimum colormap value
# If only SDE output is needed, set g4_mat = None
def plot_slice(
    g4_mat,
    sde_mat,
    axis="x",
    minval=None,
    xmin=0,
    xmax=20,
    ymin=-10,
    ymax=10,
    zmin=-10,
    zmax=10,
    voxelOrigin=np.array([0, -10, -10]),
):

    sde_mat = define_ROI(sde_mat, xmin, xmax, ymin, ymax, zmin, zmax, voxelOrigin)
    if g4_mat is not None:
        g4_mat = define_ROI(g4_mat, xmin, xmax, ymin, ymax, zmin, zmax, voxelOrigin)

    if axis == "x":
        idx = 0
        labels = ("y [cm]", "z [cm]")
        extent = [ymin, ymax, zmin, zmax]
    elif axis == "y":
        idx = 1
        labels = ("z [cm]", "x [cm]")
        extent = [xmin, xmax, zmin, zmax]
    else:
        idx = 2
        labels = ("y [cm]", "x [cm]")
        extent = [xmin, xmax, ymin, ymax]

    slicer = lambda m, i: [
        m[m.shape[0] // 2, :, :],
        m[:, m.shape[1] // 2, :],
        m[:, :, m.shape[2] // 2],
    ][i]
    sde_slice = slicer(sde_mat, idx)
    if g4_mat is not None:
        g4_slice = slicer(g4_mat, idx)

    with np.errstate(divide="ignore"):
        sde_log = np.log10(sde_slice)
        if g4_mat is not None:
            g4_log = np.log10(g4_slice)
    sde_log[np.isneginf(sde_log)] = np.nan
    if g4_mat is not None:
        g4_log[np.isneginf(g4_log)] = np.nan

    if minval is not None:
        vmin = minval
    else:
        vmin = np.nanmin([sde_log] + ([g4_log] if g4_mat is not None else []))
    vmax = np.nanmax([sde_log] + ([g4_log] if g4_mat is not None else []))

    if g4_mat is not None:
        f, ax = plt.subplots(1, 2, sharey=True, figsize=(9.2, 5), layout="compressed")
        for i, (data, title) in enumerate(zip([g4_log, sde_log], ["Geant4", "SDE"])):
            im = ax[i].imshow(
                data.T, extent=extent, origin="lower", vmin=vmin, vmax=vmax
            )
            ax[i].set_title(title, fontweight="bold")
            ax[i].set_xlabel(labels[1])
            ax[i].minorticks_on()

            if i == 0:
                ax[i].set_ylabel(labels[0])
        f.colorbar(
            im,
            ax=ax.ravel().tolist(),
            label=r"Log$_{10}$(Dose per primary) [MeV/g]",
            pad=0.02,
        )
    else:
        f, ax = plt.subplots(figsize=(7.2, 5), layout="compressed")
        im = ax.imshow(sde_log.T, extent=extent, origin="lower", vmin=vmin, vmax=vmax)
        f.colorbar(im, label=r"Log$_{10}$(Dose per primary) [MeV/g]", pad=0.02)
        ax.set_xlabel(labels[1])
        ax.set_ylabel(labels[0])
        ax.minorticks_on()

    return f


# Compares 2D projections
# set a minval to change the minimum colormap value
def plot_projection(
    g4_mat,
    sde_mat,
    minval=None,
    xmin=0,
    xmax=20,
    ymin=-10,
    ymax=10,
    zmin=-10,
    zmax=10,
    voxelOrigin=np.array([0, -10, -10]),
):
    g4_mat = define_ROI(g4_mat, xmin, xmax, ymin, ymax, zmin, zmax, voxelOrigin)
    sde_mat = define_ROI(sde_mat, xmin, xmax, ymin, ymax, zmin, zmax, voxelOrigin)
    g4_dose_proj = np.sum(g4_mat, axis=2)
    sde_dose_proj = np.sum(sde_mat, axis=2)

    # For a common colorbar ignoring infinite values when taking the log
    with np.errstate(divide="ignore"):
        g4_log = np.log10(g4_dose_proj)
        sde_log = np.log10(sde_dose_proj)

    g4_log[np.isneginf(g4_log)] = np.nan
    sde_log[np.isneginf(sde_log)] = np.nan

    if minval is not None:
        vmin = minval
    else:
        vmin = np.nanmin([g4_log, sde_log])

    vmax = np.nanmax([g4_log, sde_log])

    f, ax = plt.subplots(1, 2, sharey=True, figsize=(9.2, 5), layout="compressed")
    im1 = ax[0].imshow(
        g4_log.T, extent=[xmin, xmax, ymin, ymax], origin="lower", vmin=vmin, vmax=vmax
    )
    im2 = ax[1].imshow(
        sde_log.T, extent=[xmin, xmax, ymin, ymax], origin="lower", vmin=vmin, vmax=vmax
    )
    f.colorbar(
        im2,
        ax=ax.ravel().tolist(),
        label=r"Log$_{10}$(Integrated dose per primary) [MeV/g]",
        pad=0.02,
    )
    ax[0].set_ylabel("y [cm]")
    ax[0].set_xlabel("x [cm]")
    ax[1].set_xlabel("x [cm]")
    ax[0].set_title("Geant4", fontweight="bold")
    ax[1].set_title("SDE", fontweight="bold")
    plt.minorticks_on()
    return f


# Compares multiple 1D dose distributions either summed (projection) or from a slice
# Define if projected Bragg peak or from a slice
# Define names of the curves for other mats
# Change name of reference if not Geant4 (for example "QGSP_BIC_EMZ")
# maxdif is the upper ylim of the difference plot
def plot_multiple_bragg_peaks(
    g4_mat,
    *other_mats,
    how="proj",
    names=None,
    ref_name="Geant4",
    xmin=0,
    xmax=20,
    ymin=-10,
    ymax=10,
    zmin=-10,
    zmax=10,
    voxelOrigin=np.array([0, -10, -10]),
    maxdif=None,
):
    g4_mat = define_ROI(g4_mat, xmin, xmax, ymin, ymax, zmin, zmax, voxelOrigin)
    others = [
        define_ROI(m, xmin, xmax, ymin, ymax, zmin, zmax, voxelOrigin)
        for m in other_mats
    ]
    if names is None:
        names = [f"Mat{i + 1}" for i in range(len(others))]
    x = np.linspace(xmin, xmax, g4_mat.shape[0])
    f, (ax, ax_diff) = plt.subplots(
        2,
        1,
        figsize=(7.5, 5),
        layout="compressed",
        gridspec_kw={"height_ratios": [3, 1]},
        sharex=True,
    )
    if how == "proj":
        g4_dose1d = np.sum(g4_mat, axis=(1, 2))
    else:
        g4_dose1d = g4_mat[:, g4_mat.shape[1] // 2, g4_mat.shape[2] // 2]

    for i, (mat, name) in enumerate(zip(others, names)):
        if how == "proj":
            ax.set_ylabel(r"Total dose per primary [MeV/g]")
            dose1d = np.sum(mat, axis=(1, 2))
        else:
            ax.set_ylabel(r"Dose per primary [MeV/g]")
            dose1d = mat[:, mat.shape[1] // 2, mat.shape[2] // 2]
        ax.plot(x, dose1d, label=name)
        ax_diff.plot(
            x,
            compute_percentage_difference(g4_dose1d, dose1d),
            label=f"{name} vs {ref_name}",
        )
    ax.plot(x, g4_dose1d, "--", color="black", label=ref_name)
    ax.legend()
    ax.set_xlabel("Depth [cm]")
    ax.minorticks_on()
    if how == "proj":
        ax.ticklabel_format(style="sci", axis="y", scilimits=(3, 3))

    # Difference
    ax_diff.set_ylabel("Diff vs ref [%]")
    ax_diff.set_xlabel("Depth [cm]")
    ax_diff.minorticks_on()
    if maxdif is not None:
        ax_diff.set_ylim(0, maxdif)
    # ax_diff.legend()
    return f


# Compares lateral profiles at several depth values (cuts 1, 2 and 3 are indices in x)
# Only show profile between lowlim and uplim horizontal limits
# If only SDE output is needed, set g4_mat = None
def plot_lateral_profiles(
    g4_mat,
    sde_mat,
    xmin=0,
    xmax=20,
    ymin=-10,
    ymax=10,
    zmin=-10,
    zmax=10,
    voxelOrigin=np.array([0, -10, -10]),
    cuts=None,
    lowlim=-1.5,
    uplim=1.5,
):
    if cuts is None:
        cuts = [30, 50, 75]

    sde_mat = define_ROI(sde_mat, xmin, xmax, ymin, ymax, zmin, zmax, voxelOrigin)
    if g4_mat is not None:
        g4_mat = define_ROI(g4_mat, xmin, xmax, ymin, ymax, zmin, zmax, voxelOrigin)
        z = np.linspace(zmin, zmax, g4_mat.shape[2])
    else:
        z = np.linspace(zmin, zmax, sde_mat.shape[2])

    depths = [cut / 10 for cut in cuts]  # cm
    f, ax = plt.subplots(figsize=(7.2, 5), layout="compressed")
    used_colors = []

    for cut in cuts:
        if g4_mat is not None:
            (line,) = ax.plot(z, g4_mat[cut, g4_mat.shape[1] // 2, :])
            color = line.get_color()
            ax.scatter(
                z, sde_mat[cut, sde_mat.shape[1] // 2, :], color=color, marker="^"
            )
        else:
            (line,) = ax.plot(z, sde_mat[cut, sde_mat.shape[1] // 2, :])
            color = line.get_color()
        used_colors.append(color)

    # Depth legend (always shown)
    depth_lines = [
        Line2D([0], [0], color=c, lw=2, label=f"{d:.1f} cm")
        for c, d in zip(used_colors, depths)
    ]
    first_legend = ax.legend(handles=depth_lines, title="Depths", loc="upper right")
    ax.add_artist(first_legend)

    # Only add model legend if Geant4 is present
    if g4_mat is not None:
        model_lines = [
            Line2D([0], [0], color="k", lw=2, label="Geant4"),
            Line2D([0], [0], marker="^", color="k", linestyle="None", label="SDE"),
        ]
        ax.legend(handles=model_lines, loc="upper left")

    ax.set_xlim(lowlim, uplim)
    ax.set_xlabel("z [cm]")
    ax.set_ylabel(r"Dose per primary [MeV/g]")
    ax.minorticks_on()
    return f


# Shows all central slices for one 3D dose distribution only
def plot_all_slices(
    mat,
    xmin=0,
    xmax=20,
    ymin=-10,
    ymax=10,
    zmin=-10,
    zmax=10,
    voxelOrigin=np.array([0, -10, -10]),
):
    mat = define_ROI(mat, xmin, xmax, ymin, ymax, zmin, zmax, voxelOrigin)
    slicex = mat[mat.shape[0] // 2, :, :]
    slicey = mat[:, mat.shape[1] // 2, :]
    slicez = mat[:, :, mat.shape[2] // 2]
    # For a common colorbar ignoring infinite values when taking the log
    with np.errstate(divide="ignore"):
        slicex_log = np.log10(slicex)
        slicey_log = np.log10(slicey)
        slicez_log = np.log10(slicez)

    for s in (slicex_log, slicey_log, slicez_log):
        s[np.isneginf(s)] = np.nan

    mins = [
        np.nanmin(s)
        for s in (slicex_log, slicey_log, slicez_log)
        if np.isfinite(s).any()
    ]
    maxs = [
        np.nanmax(s)
        for s in (slicex_log, slicey_log, slicez_log)
        if np.isfinite(s).any()
    ]

    if mins and maxs:
        vmin, vmax = min(mins), max(maxs)
    else:
        vmin, vmax = 0.0, 1.0

    f, ax = plt.subplots(1, 3, figsize=(14, 5), layout="compressed")
    im1 = ax[0].imshow(
        slicex_log.T,
        extent=[ymin, ymax, zmin, zmax],
        origin="lower",
        vmin=vmin,
        vmax=vmax,
    )
    im2 = ax[1].imshow(
        slicey_log.T,
        extent=[xmin, xmax, zmin, zmax],
        origin="lower",
        vmin=vmin,
        vmax=vmax,
    )
    im3 = ax[2].imshow(
        slicez_log.T,
        extent=[xmin, xmax, ymin, ymax],
        origin="lower",
        vmin=vmin,
        vmax=vmax,
    )
    f.colorbar(
        im3,
        ax=ax.ravel().tolist(),
        label=r"Log$_{10}$(Dose per primary) [MeV/g]",
        pad=0.02,
    )
    ax[0].set_title("X central slice", fontweight="bold")
    ax[1].set_title("Y central slice", fontweight="bold")
    ax[2].set_title("Z central slice", fontweight="bold")
    for a in ax:
        a.set_aspect("equal")

    ax[0].set_xlabel("y [cm]")
    ax[0].set_ylabel("z [cm]")

    ax[1].set_xlabel("x [cm]")
    ax[1].set_ylabel("z [cm]")

    ax[2].set_xlabel("x [cm]")
    ax[2].set_ylabel("y [cm]")

    ax[0].minorticks_on()
    ax[1].minorticks_on()
    ax[2].minorticks_on()
    return f


# Performs the gamma analysis using Pymedphys package (Wendling et al 2007 http://dx.doi.org/10.1118/1.2721657)
# This function uses a fast linear interpolation, with a default fraction of 10.
# This means that the evaluation grid is interpolated at a step size of distance_mm_threshold/10.
def pymedphys_gamma(
    ref_img,
    target_img,
    dta=2,
    dd=1,
    is_local=True,
    th_percent=1,
    xmin=0,
    xmax=20,
    ymin=-10,
    ymax=10,
    zmin=-10,
    zmax=10,
    voxelOrigin=np.array([0, -10, -10]),
    plot_img=False,
    vmax=None,
):
    ref_img = define_ROI(ref_img, xmin, xmax, ymin, ymax, zmin, zmax, voxelOrigin)
    target_img = define_ROI(target_img, xmin, xmax, ymin, ymax, zmin, zmax, voxelOrigin)
    x = np.arange(0, (ymax - ymin) * 10, 1)
    y = np.arange(0, (xmax - xmin) * 10, 1)
    z = np.arange(0, (zmax - zmin) * 10, 1)
    coords = (y, x, z)
    gamma_image = pymedphys.gamma(
        coords,
        ref_img,
        coords,
        target_img,
        dose_percent_threshold=dd,
        distance_mm_threshold=dta,
        lower_percent_dose_cutoff=th_percent,
        local_gamma=is_local,
    )

    n_valid = (~np.isnan(gamma_image)).sum()
    n_passes = (gamma_image < 1).sum()
    pass_rate = 100.0 * n_passes / n_valid

    if is_local:
        method = "local"
    else:
        method = "global"

    print(
        f"{pass_rate:.2f}% pass ({method}) - DTA = {dta} mm, DD = {dd}%, "
        f"Threshold = {th_percent}%"
    )
    g = np.nan_to_num(gamma_image, nan=0)
    if plot_img:
        im = g[:, :, g.shape[2] // 2].T
        f, ax = plt.subplots(figsize=(6, 4), layout="compressed")
        plt.title(
            f"{pass_rate:.2f}% pass ({method}) \n DTA = {dta} mm, DD = {dd}%, "
            f"Threshold = {th_percent}%",
            fontweight="bold",
        )
        if vmax is not None:
            im = ax.imshow(
                im,
                extent=[xmin, xmax, ymin, ymax],
                cmap="turbo",
                interpolation="spline16",
                aspect="auto",
                vmax=vmax,
                origin="lower",
            )
        else:
            im = ax.imshow(
                im,
                extent=[xmin, xmax, ymin, ymax],
                cmap="turbo",
                interpolation="spline16",
                aspect="auto",
                origin="lower",
            )
        ax.set_xlabel("x [cm]")
        ax.set_ylabel("y [cm]")
        if xmax < 10:
            ax.set_xticks(np.arange(0, xmax, 2))
        else:
            ax.set_xticks(np.arange(0, xmax, 5))
        cbar = plt.colorbar(im, pad=0.02)
        cbar.set_label(r"Gamma index $\gamma$")
        ax.minorticks_on()
    return pass_rate, g


# Modifies a density matrix by adding a slab that starts at tmin (mm) and finishes at tmax (mm) with a density rho
# For cases with lateral heterogeneity, define a ycut that separates top to bottom. An insert will be located at the top half by default.
# If lower is set to true, the insert will be located at the bottom half.
def include_new_material(density_mat, tmin, tmax, rho, ycut=None, lower=False):
    new_mat = density_mat.copy()
    if ycut is None:
        new_mat[tmin:tmax, :, :] = rho
    else:
        if lower:
            new_mat[tmin:tmax, ycut:, :] = rho
        else:
            new_mat[tmin:tmax, :ycut, :] = rho
    return new_mat


# From depth-dose 1D distributions, compute the percentage differences found between them
def compute_percentage_difference(g4_dd, sde_dd, max_diff=50, threshold=1e-6):
    diff = np.zeros_like(g4_dd, dtype=float)
    # Mask for small g4 values because they make the diff value go crazy after bragg peak
    small_mask = g4_dd < threshold
    diff[small_mask] = max_diff
    # Mask for normal values
    normal_mask = ~small_mask
    diff[normal_mask] = (
        100 * np.abs(g4_dd[normal_mask] - sde_dd[normal_mask]) / g4_dd[normal_mask]
    )
    # Clip values above max_diff just in case
    diff = np.clip(diff, 0, max_diff)
    return diff


# From 3D dose distribution, compute percentage difference voxel-to-voxel
# Differences are normalised by percentile (default: 99th percentile of ref dose. If None, local difference is used)
def compute_percentage_voxelDiff(
    ref_img,
    target_img,
    xmin,
    xmax,
    ymin,
    ymax,
    zmin,
    zmax,
    voxelOrigin=np.array([0, -10, -10]),
    th_percent=1,
    norm_percentile=99,
):
    dose_ref = define_ROI(ref_img, xmin, xmax, ymin, ymax, zmin, zmax, voxelOrigin)
    dose_eval = define_ROI(target_img, xmin, xmax, ymin, ymax, zmin, zmax, voxelOrigin)

    mask = dose_ref > np.max(dose_ref) * th_percent / 100

    if norm_percentile is None:
        d_norm = dose_ref[mask]
    else:
        d_norm = np.percentile(dose_ref[mask], norm_percentile)

    diff = np.zeros_like(dose_ref, dtype=float)
    diff[mask] = 100 * np.abs(dose_ref[mask] - dose_eval[mask]) / d_norm

    img = diff[:, :, diff.shape[2] // 2].T
    f, ax = plt.subplots(figsize=(6, 4), layout="compressed")
    im = ax.imshow(
        img,
        extent=[xmin, xmax, ymin, ymax],
        cmap="nipy_spectral",
        interpolation="spline16",
        aspect="auto",
        origin="lower",
    )
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("y [cm]")
    if xmax < 10:
        ax.set_xticks(np.arange(0, xmax, 2))
    else:
        ax.set_xticks(np.arange(0, xmax, 5))
    cbar = plt.colorbar(im, pad=0.02)
    cbar.set_label(r"Relative difference [%]")
    ax.minorticks_on()
    return f


# Finds the proton range from an IDD profile called dose. Allows for changing from R90 to whatever percentage.
# Interpolation of the dose profile is recommended before using this function, to get a more accurate range.
def range_finder(x, dose, percentage_falloff=0.9):
    peak_idx = np.argmax(dose)
    peak_dose = dose[peak_idx]
    dose_after_peak = dose[peak_idx:]
    x_after_peak = x[peak_idx:]
    range_dose = percentage_falloff * peak_dose
    closest_idx = np.argmin(np.abs(dose_after_peak - range_dose))
    proton_range = x_after_peak[closest_idx]
    # proton_range_idx = closest_idx + peak_idx
    return proton_range  # , proton_range_idx


# Performs a comparison of proton ranges from 3D arrays, prints out raw difference
def range_comparison(
    g4_dose, sde_dose, pfalloff=0.9, total_depth=20, interp_spacing=0.01
):
    x = np.linspace(0, total_depth, total_depth * 10)
    xfine = np.linspace(0, total_depth, int(total_depth / interp_spacing))
    g4_idd = np.sum(g4_dose, axis=(1, 2))
    sde_idd = np.sum(sde_dose, axis=(1, 2))
    f1 = interp1d(x, g4_idd, kind="cubic")
    g4_fine = f1(xfine)
    f2 = interp1d(x, sde_idd, kind="cubic")
    sde_fine = f2(xfine)

    r1 = range_finder(xfine, g4_fine, percentage_falloff=pfalloff)
    r2 = range_finder(xfine, sde_fine, percentage_falloff=pfalloff)
    range_dif = np.abs(r1 - r2)
    return range_dif


# Compares 1D plots for two proton energies if second energy is available
# Let sde_dose1 be the sde matrix at energy1, g4_dose1 the geant4 matrix at energy1 and so on
def compare_bragg_peaks(
    sde_dose1,
    g4_dose1,
    energy1,
    sde_dose2=None,
    g4_dose2=None,
    energy2=None,
    max_diff1=10,
    max_diff2=25,
):
    x = np.linspace(0, 20, 200)
    g4_idd1 = np.sum(g4_dose1, axis=(1, 2))
    g4_cdd1 = g4_dose1[:, 100, 100]
    sde_idd1 = np.sum(sde_dose1, axis=(1, 2))
    sde_cdd1 = sde_dose1[:, 100, 100]

    # Integrated depth-dose curve/1DProjection
    f1, (ax1, ax1_diff) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(7.5, 5),
        gridspec_kw={"height_ratios": [3, 1]},
        layout="compressed",
    )
    ax1.plot(x, sde_idd1, label=f"SDE {energy1} MeV", color="darkblue")
    ax1.plot(
        x, g4_idd1, label=f"Geant4 {energy1} MeV", linestyle="--", color="darkturquoise"
    )
    ax1.set_ylabel(r"Total dose per primary [MeV/g]")
    ax1.ticklabel_format(style="sci", axis="y", scilimits=(3, 3))
    ax1.minorticks_on()

    # Difference
    ax1_diff.plot(
        x,
        compute_percentage_difference(g4_idd1, sde_idd1),
        label=f"{energy1} MeV",
        color="darkblue",
    )
    ax1_diff.set_ylabel("Local diff [%]")
    ax1_diff.set_xlabel("Depth [cm]")
    ax1_diff.minorticks_on()
    ax1_diff.set_ylim(0, max_diff1)

    # 1DSlice
    f2, (ax2, ax2_diff) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(7.5, 5),
        gridspec_kw={"height_ratios": [3, 1]},
        layout="compressed",
    )
    ax2.plot(x, sde_cdd1, label=f"SDE {energy1} MeV", color="darkblue")
    ax2.plot(
        x, g4_cdd1, label=f"Geant4 {energy1} MeV", linestyle="--", color="darkturquoise"
    )
    ax2.set_ylabel(r"Dose per primary [MeV/g]")
    ax2.minorticks_on()

    # Difference
    ax2_diff.plot(
        x,
        compute_percentage_difference(g4_cdd1, sde_cdd1),
        label=f"{energy1} MeV",
        color="darkblue",
    )
    ax2_diff.set_ylabel("Local diff [%]")
    ax2_diff.set_xlabel("Depth [cm]")
    ax2_diff.minorticks_on()
    ax2_diff.set_ylim(0, max_diff2)

    # Plot higher energy if available
    if energy2 is not None:
        g4_idd2 = np.sum(g4_dose2, axis=(1, 2))
        g4_cdd2 = g4_dose2[:, 100, 100]
        sde_idd2 = np.sum(sde_dose2, axis=(1, 2))
        sde_cdd2 = sde_dose2[:, 100, 100]
        ax1.plot(x, sde_idd2, label=f"SDE {energy2} MeV", color="darkred")
        ax1.plot(
            x,
            g4_idd2,
            label=f"Geant4 {energy2} MeV",
            linestyle="--",
            color="darkorange",
        )
        ax1_diff.plot(
            x,
            compute_percentage_difference(g4_idd2, sde_idd2),
            label=f"{energy2} MeV",
            color="darkred",
        )
        ax2.plot(x, sde_cdd2, label=f"SDE {energy2} MeV", color="darkred")
        ax2.plot(
            x,
            g4_cdd2,
            label=f"Geant4 {energy2} MeV",
            linestyle="--",
            color="darkorange",
        )
        ax2_diff.plot(
            x,
            compute_percentage_difference(g4_cdd2, sde_cdd2),
            label=f"{energy2} MeV",
            color="darkred",
        )

    ax1.legend()
    ax2.legend()
    ax1_diff.legend()
    ax2_diff.legend()

    return f1, f2


# 1D plots only for the SDE output
def plot_bragg_peaks_SDE(sde_dose1, energy1, sde_dose2=None, energy2=None):
    x = np.linspace(0, 20, 200)

    sde_idd1 = np.sum(sde_dose1, axis=(1, 2))
    sde_cdd1 = sde_dose1[:, 100, 100]

    f1, ax1 = plt.subplots(figsize=(7.5, 4), layout="compressed")
    ax1.plot(x, sde_idd1, label=f"{energy1} MeV", color="darkblue")
    ax1.set_ylabel(r"Total dose per primary [MeV/g]")
    ax1.set_xlabel("Depth [cm]")
    ax1.legend()
    ax1.minorticks_on()

    f2, ax2 = plt.subplots(figsize=(7.5, 4), layout="compressed")
    ax2.plot(x, sde_cdd1, label=f"{energy1} MeV", color="darkblue")
    ax2.set_ylabel(r"Dose per primary [MeV/g]")
    ax2.set_xlabel("Depth [cm]")
    ax2.legend()
    ax2.minorticks_on()

    # Plot higher energy if available
    if sde_dose2 is not None and energy2 is not None:
        sde_idd2 = np.sum(sde_dose2, axis=(1, 2))
        sde_cdd2 = sde_dose2[:, 100, 100]

        ax1.plot(x, sde_idd2, label=f"{energy2} MeV", color="darkred")
        ax2.plot(x, sde_cdd2, label=f"{energy2} MeV", color="darkred")
        ax1.legend()
        ax2.legend()

    return f1, f2
