"""
This program attempts to plot the mineralogy of a terrestrial planet
in a sensible way.
Copyright (C) 2021-2023  Fabian L. Seidler, Haiyang S. Wang

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from astropy.table import Table
import os, sys
from Utilities import Planet

path = os.path.dirname(__file__)
sys.path.append(path)

# ===== Constants ====#
Rearth = 6.371e6
Mearth = 5.974e24

# ====================#

#==== PREM density profile  =====#
tmp = Table.read(path+'/data/PREM_rhofile.txt', format='ascii')
prem_r=np.asarray(tmp['col1'])
prem_rho=np.asarray(tmp['col2'])



def fancy_mineralogy_plot(P, filename=None, samples=None, PREM=False, refRhofile=False):
    radius = np.append(P.R_core[:-1], P.R_mantle[:-1]) / Rearth
    density = np.append(P.RHO_core, P.RHO_mantle) / 1e3

    fig, ax = plt.subplots(
        1, 3, figsize=(10, 10), gridspec_kw={"width_ratios": [3, 1, 1]}, sharey=True
    )

    ax2 = plot_mineralogy(P, PREM=PREM, refRhofile=refRhofile, ax=ax[0])
    ax_pressure = ax[1]
    ax_temperature = ax[2]

    # pressure and temperature
    # --------------------------

    best_fit_pressure = np.append(P.P_core, P.P_mantle)
    if P.CMF == 0:
        best_fit_pressure = np.append(best_fit_pressure[0], best_fit_pressure)
    if P.MMF == 0:
        best_fit_pressure = np.append(best_fit_pressure, best_fit_pressure[-1])
    best_fit_radii = np.append(P.R_core, P.R_mantle)
    ax_pressure.plot(
        best_fit_pressure / 1e9,
        best_fit_radii / Rearth,
        "b",
        label="best fit",
        linewidth=3,
        zorder=1,
    )
    if samples is not None:
        for i in range(samples.shape[0]):
            tmp = samples[i]
            # try:
            ax_pressure.plot(
                tmp[1] / 1e9, tmp[0] / Rearth, "lightblue", alpha=0.3, zorder=-1
            )
            # except:
            # pass
    ax_pressure.yaxis.set_tick_params(labelleft=True)
    ax_pressure.xaxis.set_tick_params(top=True, labeltop=True)
    ax_pressure.set_xlabel("Pressure [GPa]")
    ax_pressure.legend(loc=1)
    ax_pressure.grid(True, ls="dotted", color="lightgrey")

    ax_temperature.plot(
        np.append(np.append(P.T_core[0], P.T_core), P.T_mantle),
        np.append(P.R_core, P.R_mantle) / Rearth,
        "r",
        linewidth=3,
        zorder=1,
        label="best fit",
    )
    if samples is not None:
        for i in range(samples.shape[0]):
            tmp = samples[i]
            # try:
            ax_temperature.plot(
                tmp[2], tmp[0][1:] / Rearth, "salmon", alpha=0.3, zorder=-1
            )
            # except:
            #    pass
    ax_temperature.yaxis.set_tick_params(labelleft=True)
    ax_temperature.xaxis.set_tick_params(top=True, labeltop=True)
    ax_temperature.set_xlabel("Temperature [K]")
    ax_temperature.legend(loc=1)
    ax_temperature.grid(True, ls="dotted", color="lightgrey")

    # plot CMB in pressure and temperature plot as well
    for _ax in [ax_pressure, ax_temperature]:
        _ax.axhline(
            P.R_core[-1] / Rearth, color="k", ls="dotted", linewidth=1, zorder=-2
        )

    if filename is not None:
        plt.savefig(filename + ".pdf")

    return fig, ax, ax2


def plot_mineralogy(P, PREM=False, refRhofile=False, ax=None, filename=None):
    radius = np.append(P.R_core[:-1], P.R_mantle[:-1]) / Rearth
    density = np.append(P.RHO_core, P.RHO_mantle) / 1e3
    pressure = np.array(P.mantle_dict["pressure"]) * 1e5
    mineralogy = P.mineralogy

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 11))
    ax1 = ax
    ax2 = ax1.twiny()

    #### plot density
    ax2.plot(density, radius, color="k", linewidth=3, zorder=1)

    ####plot prem density
    if PREM==True:
        ax2.plot(prem_rho, prem_r+(max(radius)-max(prem_r)), color='r', linewidth=2, zorder=1)
    if refRhofile==True:
        ### plot modelled Earth as reference
        mantle_chemsys = {"SiO2": 0.45,
            "MgO": 0.378,
            "FeO": 0.0805,
            "Al2O3": 0.0445,
            "CaO": 0.0355,
            "Na2O": 0.0036}   ## MS1955's mantle composition (which will be automatically normalised to 1 in the interior_model)

        core_chemsys = {"Fe": 82.8, "S": 1.83, "Ni": 5.06, "Si": 4.96} ## Wang et al. 2018 core composition (which will be automatically normalised to 1 in the interior_model)
        R = 1
        CMF = 0.325

        PE = Planet(R, CMF, mantle_chemsys, core_chemsys, P_max=140e9)
        PE.solve_structure()
        Eradius = np.append(PE.R_core[:-1], PE.R_mantle[:-1]) / Rearth
        Edensity = np.append(PE.RHO_core, PE.RHO_mantle)/1e3
        ax2.plot(Edensity, Eradius+(max(radius)-max(Eradius)), color='k', linewidth=2, linestyle='dashed', zorder=1)


    #### plot core
    _tmp = np.zeros_like(radius)
    _tmp[: P.N_mantle + 1] = 100
    ax1.fill_betweenx(
        radius, np.zeros(P.N_core + P.N_mantle), _tmp, color="grey", zorder=-1
    )
    _x = (np.amin(radius[P.N_core :])) / 2
    core_text = "Fe-Ni-Si-S core"
    ax1.annotate(core_text, (40, _x), zorder=1)

    # xlim and ylim
    ax1.set_ylim([0, np.amax(radius)])
    ax1.set_xlim([0, 100])
    ax2.set_ylim([0, np.amax(radius)])
    ax2.set_xlim([np.amin(density)-1, np.amax(density)+1]) #9 for 0.5R planets; 13 for 1R planets; #20 - for 1.5R planets

    # labels
    mode = "volume"
    if mode == "volume" or mode == "Volume" or mode == "vol" or mode == "Vol":
        mode_str = "[vol]"
    elif mode == "molar" or mode == "mol":
        mode_str = "[mol]"
    ax1.set_ylabel(r"Radius [$R_\oplus$]")
    ax1.set_xlabel(r"Phase % " + mode_str)
    ax1.xaxis.set_label_position("top")
    ax1.xaxis.tick_top()
    ax2.set_xlabel(r"Density [$\mathrm{g/cm}^3$]")
    ax2.xaxis.set_label_position("bottom")
    ax2.xaxis.tick_bottom()

    # ticks
    ax2.minorticks_on()
    ax1.minorticks_on()

    # define colors
    colors = {
        "Wus": "sandybrown",
        "Sp": "maroon",
        "P": "tomato",
        "Pp": "r",
        "Wad": "greenyellow",
        "O": "darkkhaki",
        "Ring": "yellowgreen",
        "Aki": "plum",
        "Pl": "lightcoral",
        "st": "lightpink",
        "Cpx": "aliceblue",
        "Opx": "lightsteelblue",
        "CF": "orange",
        "C2/c": "lavender",
        "ca-p": "khaki",
        "ky": "b",
        "coe": "aquamarine",
        "seif": "magenta",
        "Gt": "moccasin",
        "q": "paleturquoise",
        "neph": "tab:blue",
    }
    mineral_names = {
        "Wus": "wus",
        "P": "mg-pv",
        "Pp": "mg-postpv",
        "Wad": "wad",
        "O": "ol",
        "Ring": "ring",
        "Aki": "aki",
        "Pl": "plg",
        "st": "stv",
        "Cpx": "cpx",
        "Opx": "opx",
        "CF": "cf",
        "C2/c": "hp-cpx",
        "ky": "ky",
        "ca-p": "ca-pv",
        "coe": "coe",
        "seif": "seif",
        "Gt": "gt",
        "q": "q",
        "Sp": "sp",
        "neph": "neph",
    }
    bottom = np.zeros(P.N_mantle + 1)

    #### plot minerals/mantle
    for i in range(len(P.mineral_phases)):
        mineralfit = interp1d(
            pressure,
            mineralogy[:, i],
            kind="slinear",
            bounds_error=False,
            fill_value="extrapolate",
        )
        mineral = mineralfit(P.P_mantle)
        RADIUS = P.R_mantle / Rearth
        ax1.fill_betweenx(
            RADIUS,
            bottom,
            bottom + mineral,
            step="mid",
            label=P.mineral_phases[i],
            edgecolor=None,
            facecolor=colors[P.mineral_phases[i]],
        )

        #### labeling the minerals
        tmp = mineral.copy()
        tmp[np.where(mineral != 0)] = 1
        tmp2 = tmp * RADIUS
        _x = np.mean(tmp2[np.where(tmp2 != 0)])  # + radius[P.N_core]
        tmp_index = np.argmin(abs(RADIUS - _x))
        _y = bottom[tmp_index] + (mineral[tmp_index]) / 2
        if mineral[tmp_index] <= 15:
            rotation = -90
        else:
            rotation = 0
        volume = np.sum(mineral[:-1] * np.diff(RADIUS))
        if volume < 0.5:
            if _y < 40:
                _y_arrow = _y + 10 * np.random.rand()
            else:
                _y_arrow = _y - 10 * np.random.rand()
            ax1.annotate(
                mineral_names[P.mineral_phases[i]],
                xy=(_y, _x),
                xytext=(_y_arrow, _x - 0.03),
                arrowprops=dict(facecolor="black", arrowstyle="-"),
            )
        else:
            ax1.annotate(
                mineral_names[P.mineral_phases[i]], (_y, _x), rotation=rotation
            )

        ####
        bottom += mineral

    # plt.tight_layout()

    if filename is not None:
        plt.savefig(filename + ".pdf")

    return ax2
