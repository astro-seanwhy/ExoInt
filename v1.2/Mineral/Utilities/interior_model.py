"""
This program evaluates the interior structure of a terrestrial planet of given radius,
based on the mineralogy calculated by Perple_X and the liquid iron equation
of state from Kuwayama et al. 2020.
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
from scipy.optimize import brentq
from scipy.interpolate import UnivariateSpline, interp1d
import sys, os
import random, time
from run_perplex import calculate_mineralogy, make_mineralogy_grid
import warnings
from math import pow, exp

path = os.path.dirname(__file__)

# ====== Constants ========#
kB = 1.38e-23  # J K⁻¹
u = 1.6605e-27  # kg
Avogadro = 6.022e23
R = 8.31446  # m³ Pa K⁻¹ mol⁻¹



G = 6.674e-11
Mearth = 5.974e24
Rearth = 6.371e6

mO = 15.994
mMg = 24.305
mSi = 28.0855
mFe = 55.847
mAl = 26.981
mCa = 40.078
mNa = 22.9898

mSiO2 = mSi + 2 * mO
mMgO = mMg + mO
mFeO = mFe + mO
mAl2O3 = 2 * mAl + 3 * mO
mCaO = mCa + mO
mNa2O = 2 * mNa + mO

atomic_mass = {
    "Fe": 55.845,
    "Ni": 58.6934,
    "Si": 28.0867,
    "O": 15.994,
    "S": 32.065,
    "C": 12.011,
    "H": 1.00794,
    "U": 238.0289,
    "Th": 232.0381,
    "Cr": 51.9951,
    "Mn": 54.93805,
    "P": 30.973762,
    "Co": 58.9332,
    "Na": 22.989769,
    "Mg": 24.305,
    "Ca": 40.078,
    "Al": 26.981539,
}


Fe_atomic_mass = atomic_mass['Fe'] * u
Fe_molar_mass = Fe_atomic_mass * Avogadro


# ======= planetary structure equations ========#
def update_radius(m, R, RHO):
    """
    For a planet dissected into nested shells, update the radii of
    each shell given new densities RHO.

    Parameters
    ----------
        m : np.ndarray
            array with the (fixed) masses of each subshell
        R : np.array
            radii of each subshell
        RHO : np.ndarray
            array with the densities of each subshell
    Returns
    -------
        r : np.ndarray
            new radii of each subshell (same shape as R)
    """
    r = R.copy()
    for i in range(len(RHO)):
        r[i + 1] = np.cbrt(
            3 / (4 * np.pi) * (m[i] / RHO[i] + 4 * np.pi / 3 * r[i] ** 3)
        )
    return r


def update_pressure(M, R, RHO, Psurf, G=6.676e-11):
    """
    Evaluates the hydrostatic equation.

    The integration method is Explicit Euler. Hence, this function is valid if M_i ~ M_i+1, R_i ~ R_i+1,
    or equivalently if the masses of shells m_i are sufficiently small.

    Parameters
    ----------
        M : np.ndarray
            cumulative masses, in kg
        R : np.ndarray 
            radii of each subshell, in m
        RHO : np.ndarray
           densities of each subshell, in kg/m^3
        Psurf : float
            surface pressure
    Returns
    -------
        P : np.ndarray
            pressures, same shape as R
    """
    P = np.ones_like(M) * Psurf
    for i in range(len(P) - 2, -1, -1):
        if R[i] > 0:
            P[i] = P[i + 1] + G * M[i] / (R[i] ** 2) * RHO[i - 1] * (R[i + 1] - R[i])
    if R[0] == 0.0:
        P[0] = P[1] + 2 / 3 * np.pi * G * RHO[0] ** 2 * R[1] ** 2
    return P


def update_temperature(P, KS, gamma, T_surface):
    """
    Computes an adiabatic gradient along a pressure array.

    Parameters
    ----------
        P : np.ndarray
            pressures, unit Pa
        KS : np.ndarray
            adiabatic bulk modulus, unit Pa
        gamma : np.ndarray
            Grüneisen-parameter in each layer, dimensionless
        T_surface : float
            surface temperature of layer, unit K
    Returns
    -------
        T : np.ndarray
            Temperature, same shape as KS & gamma.
    """

    T = np.ones(len(P) - 1) * T_surface
    for i in range(len(T) - 2, -1, -1):
        T[i] = T[i + 1] - gamma[i] / KS[i] * T[i] * (P[i + 1] - P[i])
    return T


def geothermal_gradient(P, T_surf):
    """
    Estimates an adiabatic mantle temperature profile based on a simulation
    with the burnman package (), where I used the mineralogy of Earth and the
    adiabatic equation and started to integrate along increasing pressures 
    starting from various surface temperatures, after which I applied a simple fit.
    It should be accurate for planets of < 2 Mearth.

    Parameters
    ----------
        P : float or np.ndarray
            pressure(s) of interest, in Pa
        T_surf : float
            surface temperature
    Returns
    -------
        T : float or np.dnarray
            Temperatures, in K
    """
    a = 7.07890407e-12 * T_surf**2 - 1.26873023e-08 * T_surf + 7.55756985e-06
    n = 0.7607342014182283 - 3.5920125989830075e-05 * (T_surf - 1000)
    return a * P**n + T_surf




def generate_pressure_grid(P_surface, P_max, dP=[0.75e9, 2e9]):
    """
    This function is intended to create a pressure profile along which the thermodynamics 
    will be evaluated with Gibbs-free-energy-minimization (GFEM). Since each step with GFEM 
    takes some runtime, we want to do it in as little steps as possible while still maintaining 
    sufficient accuracy. Hence, this function is invoked. In its simplest form, as presented here, 
    it will just produce an upper mantle with more resolution (since there more mineral transitions happen),
    and a lower mantle with lower resolution.

    Parameters
    ----------
        P_surface : float
            surface pressure, in Pa
        P_max : float
            maximum pressure; estimated in the function, in Pa
        dP : list
            resultion of the mantle layers; [dP_upper_mantle, dP_lower_mantle], in Pa
    Returns
    -------
        P : np.ndarray
            pressure array, in Pa
    """

    if P_surface < 30e9:
        if P_max > 30e9:
            P_upper_mantle = np.linspace(
                P_surface, 30e9, int((30e9 - P_surface) / dP[0])
            )
            P_lower_mantle = np.linspace(30e9, P_max, int((P_max - 30e9) / dP[1]))

            P = np.hstack((P_upper_mantle, P_lower_mantle[1:]))
        else:
            N = max(2, int((P_max - P_surface) / dP[0]))
            P = np.linspace(P_surface, P_max, N)

    else:
        N = max(2, int((P_max - P_surface) / dP[1]))
        P = np.linspace(P_surface, P_max, N)

    return P



def call_mineralogy_calculation(
    P_surface, P_CMB, T_surface, mantle_chemsys, subpath_name=None
):
    """
    Calculates the mineralogy with Perple_X

    Parameters
    ----------
        P_surface : float
            surface pressure, in Pa
        P_CMB : float
            pressure at core-mantle boundary, in Pa
        T_surface : float
            surface temperature, in K
        mantle_chemsys : dict
            dictionary containing the abundances (in massfractions) of SiO2, MgO, FeO, Al2O3, CaO, Na2O
    Returns
    -------
        mantle_dict : dict
            dictionary containing mantle properties:
                {'temperature','pressure','density','Grueneisen','KS', 'speeds', 'phases'}
                Each of them is a a numpy array, containing the respective parameters.
                'speeds' is a 2d-array and contains the seismic speeds
                'phases' is a 2d-array with volume abundance of mineral phases as entries;
                the columns correspond to the respective phase listed in 'Phases' (see below).
        Phases: list
                Names of the mineral phases.

    """
    if subpath_name is None:
        random.seed(mantle_chemsys.get("SiO2"))
        subpath_name = str(
            random.randint(a=10000, b=99999)
        )  # this is needed to run PerPleX
        while os.path.isdir(path + "/tmp_files/" + subpath_name):
            subpath_name = str(random.randint(a=10000, b=99999))
    os.system("mkdir " + path + "/tmp_files/" + subpath_name)


    # generate the pressure and temperature grid
    P = generate_pressure_grid(P_surface, P_CMB)
    T = geothermal_gradient(P, T_surface)

    mantlegrid = open(path + "/tmp_files/" + subpath_name + "/mantlegrid.dat", "w")
    for i in range(len(P)):
        mantlegrid.write(str(P[i] / 1e5) + " " + str(T[i]) + "\n")
    mantlegrid.close()

    calculate_mineralogy(mantle_chemsys, path + "/tmp_files/" + subpath_name)

    # check if Perple_X has completed successfully
    if not os.path.isfile(path + "/tmp_files/" + subpath_name + "/mantlephys_1.tab"):
        os.system("rm -rf " + path + "/tmp_files/" + subpath_name)
        raise PerpleX_Error(
            "Perple_X can not find a stable mineralogy for the given input. Exiting"
        )

    mantle_dict, Phases = make_mineralogy_grid(path + "/tmp_files/" + subpath_name)
    os.system(
        "mv "
        + path
        + "/tmp_files/"
        + subpath_name
        + "/mantlephys_1.tab "
        + path
        + "/mineralogy.tab"
    )
    os.system("rm -r " + path + "/tmp_files/" + subpath_name)

    return mantle_dict, Phases


# ---------- functions for core temperature (Noack & Lasbleis 2020 )----------#
def T_CMB_warm(P_CMB, FeM, X_M0=0.11):
    """
    Temperature estimate at core-mantle boundary
    (Noack & Lasbleis 2020, A&A 638, A19, Eq. 21)
    """
    return 5400 * (P_CMB / 140e9) ** 0.48 / (1 - np.log(1 - X_M0 - FeM))


def T_CMB_hot(P_CMB, FeM):
    """
    Temperature estimate at core-mantle boundary
    (Noack & Lasbleis 2020, A&A 638, A19, Eq. 20)
    """
    return 5400 * (P_CMB / 140e9) ** 0.48 / (1 - np.log(1 - FeM))



def get_FeM(wtFeO, wtMgO):
    """
    Fe-Number, following Noack & Lasbleis 2020, A&A 638, A19
    """
    NFeO = wtFeO / mFeO
    NMgO = wtMgO / mMgO
    return NFeO / NMgO


# ================== Main class ====================#
class M_Planet:
    def __init__(self, M, CMF, mantle_chemsys, core_chemsys):
        self.N_core = 500
        self.N_mantle = 500
        self.tolerance = 0.1

        # surface boundary conditions
        self.P_surface = 1e5
        self.T_surface = 1700  # 1673.15

        # which EoS should be used for the core modelling?
        self.core_eos = "Kuwayama"  ##The other option is "Seager2007"

        #print("M: ", M)
        #print("CMF: ", CMF) 

        self.M = M * Mearth
        self.CMF = CMF
        self.MMF = 1 - CMF

        self.mantle_chemsys = mantle_chemsys
        self.core_chemsys = core_chemsys

        self._check_input()

        # calculate lower core density due to light elements;
        # we do so by lowering the mean molecular weight of the material
        if self.CMF != 0.0:
            self.molar_mass_corematerial = 1 / sum(
                self.core_chemsys[core_element] / atomic_mass[core_element]
                for core_element in self.core_chemsys
            )

        self._Delta_T_CMB = 0.0  # 1400 * M**(3/4)
        self.Delta_T_CMB = None

    def _check_input(self, verbose=False):
        """
        Verifies that input parameters are correct.
        """

        if self.M <= 0:
            raise ValueError("A planet with mass<=0? You must be kidding, right?")

        if self.CMF < 0 or self.CMF > 1:
            raise ValueError("CMF<0 or CMF>1 makes no sense.")

        for entry in self.mantle_chemsys:
            if np.isnan(self.mantle_chemsys[entry]):
                self.mantle_chemsys[entry] = 0.0

        for entry in self.core_chemsys:
            if np.isnan(self.core_chemsys[entry]):
                self.core_chemsys[entry] = 0.0

        if np.any(np.asarray(list(self.mantle_chemsys.values())) < 0):
            raise ValueError(
                "Can't handle the negative abundance in mantle chemistry array. ¯\_(ツ)_/¯"
            )

        total_mantle_mass = sum(x for x in self.mantle_chemsys.values())
        if abs(1 - total_mantle_mass) > 0.001:
            if verbose:
                warnings.warn(
                    "Warning: the mass-fractions of mantle chemistry does not sum to 1! Calculation is continued - but the mantle chemistry of the planet model has been normalised to 1!"
                )

            for compound in self.mantle_chemsys:
                self.mantle_chemsys[compound] /= total_mantle_mass

        if self.CMF != 0:
            if np.any(np.asarray(list(self.core_chemsys.values())) < 0):
                raise ValueError(
                    "Don't feed me negative abundances of core light elements! Really, please don't."
                )

            if "Fe" not in self.core_chemsys:
                if np.sum(list(self.core_chemsys.values())) > 1:
                    raise ValueError(
                        "The massfractions of core elements sum to more than 1. This is unphysical, please check the input."
                    )
                else:
                    self.core_chemsys["Fe"] = 1 - sum(
                        x for x in self.core_chemsys.values()
                    )

            total_core_mass = sum(x for x in self.core_chemsys.values())
            if abs(1 - total_core_mass) > 0.001:
                if verbose:
                    warnings.warn(
                        "Warning: the mass-fractions of core chemistry do not sum to 1! They are now renormalized."
                    )

                for compound in self.core_chemsys:
                    self.core_chemsys[compound] /= total_core_mass

    def init_planet(self):
        """
        Initializes the planet as layered nested and spherical shells, each of which has
        the zero-pressure density of the respective material (mantle or core).
        """
        core_mass = self.M * self.CMF
        mantle_mass = self.M * self.MMF
        density = {"mantle": 4100.0, "core": 8300.0}
        Rc = np.cbrt(3 * core_mass / (4 * np.pi * density["core"]))
        V_mantle = mantle_mass / density["mantle"]
        R_mantle = np.cbrt(3 * V_mantle / (4 * np.pi) + Rc**3)

        # shells
        self.R_core = np.linspace(0, Rc, self.N_core + 1)
        self.R_mantle = np.linspace(Rc, R_mantle, self.N_mantle + 1)
        self.RHO_core = np.ones(self.N_core) * density["core"]
        self.RHO_mantle = np.ones(self.N_mantle) * density["mantle"]
        self.m_core = (
            4
            * np.pi
            / 3
            * (self.R_core[1:] ** 3 - self.R_core[:-1] ** 3)
            * self.RHO_core
        )
        self.m_mantle = (
            4
            * np.pi
            / 3
            * (self.R_mantle[1:] ** 3 - self.R_mantle[:-1] ** 3)
            * self.RHO_mantle
        )

        # temperature
        self.T_mantle = np.ones_like(self.RHO_mantle) * self.T_surface
        self.T_core = np.ones_like(self.RHO_core) * self.T_surface

        # adjust arrays in case CMF or MMF is 0:
        if self.CMF == 0.0:
            self.m_core = np.zeros_like(self.RHO_core)
            self.RHO_core = np.zeros_like(self.RHO_core) * np.nan
            self.P_core = np.zeros(self.N_core + 1) * np.nan
            self.T_core = np.zeros(self.N_core) * np.nan
        if self.MMF == 0.0:
            self.m_mantle = np.zeros_like(self.RHO_mantle)
            self.RHO_mantle = np.zeros_like(self.RHO_mantle) * np.nan
            self.P_mantle = np.ones(self.N_mantle + 1) * self.P_surface
            self.T_mantle = np.ones(self.N_mantle) * self.T_surface

        # calculate the cumulative mass M:
        self.M_core = np.zeros_like(self.R_core)
        if self.CMF != 0:
            for i in range(1, self.N_core + 1):
                self.M_core[i] = np.sum(self.m_core[:i])

        self.M_mantle = np.zeros_like(self.R_mantle)
        if self.MMF != 0:
            for i in range(1, self.N_mantle + 1):
                self.M_mantle[i] = np.sum(self.m_mantle[:i])
        self.M_mantle += self.M_core[-1]

        # init pressure
        if self.MMF != 0:
            self.P_mantle = update_pressure(
                self.M_mantle, self.R_mantle, self.RHO_mantle, self.P_surface
            )

        if self.CMF != 0:
            self.P_core = update_pressure(
                self.M_core, self.R_core, self.RHO_core, self.P_mantle[0]
            )

    def Msolve_structure(self, mantle_dict=None, subpath_name=None):
        """
        Solve for the interior structure and mineralogy of a planet of given mass.
        """
        self.init_planet()

        if self.successful == True and mantle_dict is not None:
            converged = False
            self.Nsteps = 0
            while converged == False:
                stored_step = self.store_current_step()
                self.step(
                    mantle_eos="perpleX",
                    core_eos=self.core_eos,
                    mantle_dict=mantle_dict,
                )
                self.Nsteps += 1
                if self.convergence_criterion(stored_step) == True:
                    converged = True
                if self.Nsteps > 100:
                    self.successful = False
                    self.error_flag = "Endless iteration."
                    converged = True
            ####

        self.CRF = self.R_core[-1] / self.R
        return self.R_mantle[-1]

    def store_current_step(self):
        """
        Stores the current step - used later to determine convergence.
        """
        if self.CMF == 0:
            RHO = self.RHO_mantle
            R = self.R_mantle
            P = self.P_mantle
            T = self.T_mantle
        elif self.MMF == 0:
            RHO = self.RHO_core
            R = self.R_core
            P = self.P_core
            T = self.T_core
        else:
            RHO = np.append(self.RHO_core, self.RHO_mantle)
            R = np.append(self.R_core, self.R_mantle)
            P = np.append(self.P_core, self.P_mantle)
            T = np.append(self.T_core, self.T_mantle)
        return {"rho": RHO, "R": R, "P": P, "T": T}

    def convergence_criterion(self, stored_step):
        """
        Determines if the calculation has converged to the specified tolerance.
        """
        if self.CMF == 0:
            RHO = self.RHO_mantle
        elif self.MMF == 0:
            RHO = self.RHO_core
        else:
            RHO = np.append(self.RHO_core, self.RHO_mantle)

        distance = abs(stored_step["rho"] - RHO)

        if np.all(distance < self.tolerance):
            return True
        else:
            return False

    def step(
        self, mantle_eos="Seager2007", core_eos="Seager2007", mantle_dict=None
    ):
        """
        Compresses the planet by one iteration.
        """
        if self.MMF != 0:
            self.P_mantle = update_pressure(
                self.M_mantle, self.R_mantle, self.RHO_mantle, self.P_surface
            )
        if self.CMF != 0:
            self.P_core = update_pressure(
                self.M_core, self.R_core, self.RHO_core, self.P_mantle[0]
            )
        else:
            self.P_core = np.ones(self.N_core) * self.P_mantle[0]

        if self.MMF != 0:
            self.update_mantle_phys(eos=mantle_eos, mantle_dict=mantle_dict)
        if self.CMF != 0:
            self.update_core_phys(core_eos)

        if self.MMF != 0:
            self.T_mantle = geothermal_gradient(self.P_mantle, self.T_surface)
        if self.CMF != 0:
            self.T_core = update_temperature(
                self.P_core,
                self.KS_core,
                self.gr_core,
                self.T_mantle[0] + self.Delta_T_CMB,
            )

        if self.CMF != 0:
            self.R_core = update_radius(self.m_core, self.R_core, self.RHO_core)
        if self.MMF != 0:
            self.R_mantle = update_radius(
                self.m_mantle, self.R_mantle, self.RHO_mantle
            ) - (self.R_mantle[0] - self.R_core[-1])
        else:
            self.R_mantle = np.ones(self.N_mantle) * self.R_core[-1]


    # calculate the thermodynamics of the mantle
    # -------------------------------------------
    def update_mantle_phys(self, eos="Seager2007", mantle_dict=None):
        """
        Computes and interpolates the profiles for the density (rho), 
        Grueneisen-parameter (gr) and adiabatic bulk modulus (KS) in the
        silicate mantle. Resulting arrays will be processed in the compression
        step (self.step)
        """
        if eos == "Seager2007":
            self.RHO_mantle = (
                4100
                + 0.00161 * (self.P_mantle[:-1] + np.diff(self.P_mantle) / 2) ** 0.541
            )
            self.gr_mantle = 1 * np.ones(self.N_mantle)
            self.KS_mantle = 1e11 * np.ones(self.N_mantle)

        elif eos == "perpleX":
            if mantle_dict is not None:
                P_thermo = np.array(mantle_dict["pressure"]) * 1e5
                RHO_mantle = np.array(mantle_dict["density"])
                gr_mantle = np.array(mantle_dict["Grueneisen"])
                KS_mantle = np.array(mantle_dict["KS"]) * 1e5

                # interpolate the pressure grid from the thermodynamic modelling to get the mantle thermodynamic properties of the planet:
                rhofunc = interp1d(P_thermo, RHO_mantle, fill_value="extrapolate")
                grfunc = interp1d(P_thermo, gr_mantle, fill_value="extrapolate")
                KSfunc = interp1d(P_thermo, KS_mantle, fill_value="extrapolate")

                self.RHO_mantle = rhofunc(
                    self.P_mantle[1:] - np.diff(self.P_mantle) / 2
                )
                self.gr_mantle = grfunc(self.P_mantle)
                self.KS_mantle = KSfunc(self.P_mantle)

    # evaluate the core equation of state
    # -----------------------------------------
    def update_core_phys(self, eos="Seager2007"):
        """
        Computes and interpolates the profiles for the density (rho), 
        Grueneisen-parameter (gr) and adiabatic bulk modulus (KS) in the
        silicate mantle. Resulting arrays will be processed in the compression
        step (self.step)
        """

        if eos == "Seager2007":
            self.RHO_core = (
                8300 + 0.00349 * (self.P_core[:-1] + np.diff(self.P_core) / 2) ** 0.528
            )

        elif eos == "Kuwayama":
            # ---------- temperature jump at the core-mantle boundary
            if type(self._Delta_T_CMB) == str:
                # Noack and Lasbleis 2020
                if self.M < 0.8*Mearth or self.M > 2*Mearth:
                    self._Delta_T_CMB = "stx14"
                    #print("M=" + str(round(self.M/Mearth, 3)) + " is out of the applicable range [0.8, 2]M_Earth of NL20!")
                    #print("Thus, the calculation of Delta_T_CMB is implemented by the stx14 approach.")

                if self._Delta_T_CMB == "NL20":
                    self.FeM = get_FeM(
                        self.mantle_chemsys["FeO"], self.mantle_chemsys["MgO"]
                    )
                    T_CMB = T_CMB_warm(
                        self.P_mantle[0], self.FeM
                    )  # get_FeM(self.mantle_chemsys['FeO'], self.mantle_chemsys['MgO']))
                    self.Delta_T_CMB = T_CMB - self.T_mantle[0]

                # Stixrude 2014
                elif self._Delta_T_CMB == "stx14":
                    self.Delta_T_CMB = 1400 * (self.M / Mearth) ** (3 / 4)
                    T_CMB = self.T_mantle[0] + self.Delta_T_CMB
                else:
                    raise NotImplementedError(
                        "No core temperature model with name '"
                        + self.Delta_T_CMB
                        + "' is implemented"
                    )

            # constant as temperature jump
            elif type(self._Delta_T_CMB) == float or type(self._Delta_T_CMB) == int:
                self.Delta_T_CMB = self._Delta_T_CMB
                T_CMB = self.T_mantle[0] + self.Delta_T_CMB
                # print(self.Delta_T_CMB)
            else:
                warnings.warn(
                    "No valid T_cmb jump method is specified; zero value of Delta_T_CMB is applied"
                )
                self.Delta_T_CMB = 0
                T_CMB = self.T_mantle[0] + self.Delta_T_CMB

            # ----------- Compute core EOS -----------#
            P0 = 1e5
            T0 = 1811.0
            V0 = 0.00000794
            KT0 = 82.1e9
            KT0dash = 5.8
            gamma0 = 2.02
            q = 0.63
            e0 = 0.68e-4
            g = -1.0
            n = 1
            rho0 = self.molar_mass_corematerial / V0 / 1e3

            def thermal_energy(x, T):
                return 3*n*R*(T + e0 * pow(x, g) * pow(T,2))

            def Delta_thermal_Pressure(x, T, rho0):
                gamma = gamma0 * pow(x,q)
                V = Fe_molar_mass*x/rho0
                return gamma/V * (thermal_energy(x,T) - thermal_energy(x,T0))

            def Cold_Pressure(x):
                return 3.*KT0 * pow(x, -2/3.) * (1-pow(x,1/3.)) * exp( 3/2.*(KT0dash-1) * (1-pow(x,1/3.)) )

            def Kuwayama_pressure(rho, T, rho0):
                """
                returns the pressure as a function of rho,T
                """
                x = rho0/rho
                P = Cold_Pressure(x) + Delta_thermal_Pressure(x,T,rho0)
                return P

            class test_params:
                def __init__(self, P, T, rho0):
                    self.P = P
                    self.T = T
                    self.rho0 = rho0 

            def f(rho, args):
                #myargs = args
                return Kuwayama_pressure(rho, args.T, args.rho0) - args.P

            def density(P, T, rho0):   
                args = test_params(P, T, rho0)
             
                rtol = 2e-15
                xtol = 2e-12
                mitr = 100
             
                xa = .99*rho0
                xb = 1e4*rho0
             
                return brentq(f, xa, xb, args, xtol, rtol, mitr)
             
            def thermal_alpha(P, T, rho0):  #expansion
                rho = density(P, T, rho0)
                alpha = -1/rho * (density(P,T+.5,rho0) - density(P,T-.5,rho0))
                return alpha 


            def alpha_from_rho(rho, T, rho0): #bulk modulus
                #rhoc = rho(P, T, rho0)
                return rho * (Kuwayama_pressure(rho+.5,T,rho0) - Kuwayama_pressure(rho-.5,T,rho0))


            #def Cv(x, T):  #isochoric_heatcapacity
            #    return 3 * n * R * (1 + 2 * e0 * x**g * T)

            # alpha is the thermal expansion coefficient
            alpha = np.zeros(self.N_core)
            KT = np.zeros(self.N_core)

            for i in range(self.N_core):
                P_inbetween = (
                    self.P_core[i] + (self.P_core[i + 1] - self.P_core[i]) / 2.0
                )
                self.RHO_core[i] = density(P_inbetween, self.T_core[i], rho0)
                alpha[i] = thermal_alpha(P_inbetween, self.T_core[i], rho0)
                
                KT[i] = alpha_from_rho(self.RHO_core[i], self.T_core[i], rho0)
            
            
            #self.gr_core = gamma0 * (self.RHO_core / rho0) ** q  ##this is wrong
            self.gr_core = gamma0 * (rho0 / self.RHO_core) ** q
            self.KS_core = (1 + self.gr_core * alpha * self.T_core) * KT


class Planet(M_Planet):
    def __init__(self, R, CMF, mantle_chemsys, core_chemsys, P_max=500e9):
        self.R = R * Rearth
        self.P_max = P_max

        super().__init__(1., CMF, mantle_chemsys, core_chemsys)

    def __foo(self, M):
        self.M = M * Mearth
        _R = self.Msolve_structure(mantle_dict=self.mantle_dict)
        return _R - self.R

    def solve_structure(self, subpath_name=None, M_low=0.001, M_max=10):
        """
        Solve for the interior structure and mineralogy of a planet of given radius.
        """

        #### call Perple_X
        self.successful = False
        self.error_flag = "Mantle thermodynamics not yet evaluated"
        if self.CMF != 1:

            # ---- compute the mantle mineralogy ----#
            start = time.time()
            mantle_dict, Phases = call_mineralogy_calculation(
                self.P_surface,
                self.P_max,
                self.T_surface,
                self.mantle_chemsys,
                subpath_name,
            )
            end = time.time()
            print("Mineralogy computation: ", end - start)
            self.mantle_dict = mantle_dict
            self.successful = True
            self.error_flag = "No error occured."

            self.mineralogy = np.asarray(self.mantle_dict["phases"])
            self.mineralogy[np.isnan(self.mineralogy)] = 0
            self.mineral_phases = Phases

        else:
            mantle_dict = None
            self.successful = True
            self.error_flag = "No error occured, but the planet has no mantle."

        if self.successful == True:
            self.M = brentq(self.__foo, M_low, M_max)
            self.P_CMB = self.P_mantle[0]
            if self.M < 0.8 or self.M > 2:
                print("Please note that the default approach -- NL20 -- of calculating Delta_T_CMB has been replaced by stx14, since the mass of " + str(round(self.M, 3)) + " has been out of the parameterised range [0.8, 1.2]M_Earth of NL20!")
