import numpy as np
import sys, os
import pandas
from astropy.table import Table
from tqdm import tqdm

from chemsysmodel import chemsysmodel
from utilities import load_data, produce_planetabuN, generate_asymmetric_random, single_simulation_report


# ======= Constants =========#
tmp = Table.read("data/atomwttc_new.txt", format="ascii")
elemid = np.asarray(tmp["elem"])


# ======== Function ==========#
def exomodel(
    stellar_abundances_file,
    Nsim=1e4,
    refsolar="A21",
    asymerr=False,
    outdir='output/',
    verbose=True,
):

    # =================== PREPARATION ==================#
    Nsim = int(Nsim)

    # has outdir the correct format?
    if outdir[-1] != '/':
        outdir += '/'

    # does the outdir directory exist?
    if not os.path.isdir(outdir):
        os.system("mkdir -p "+str(outdir))

    # specify filename (i.e., what will the output be called?)
    filename = stellar_abundances_file.split("abu_")[-1]
    filename = filename.split(".")[0]

    # load stellar abundance data
    stardex, stardexerr, stardexerrdn = load_data(
        stellar_abundances_file, asymerr=asymerr
    )

    # compute planetary abundances (as number of atoms)
    planetabuN, planetabuNerrup, planetabuNerrdn = produce_planetabuN(
        stardex, stardexerr, stardexerrdn, refsolar=refsolar
    )

    #breakpoint()

    # ========= MEAN STELLAR ABUNDANCES RESULT =============#
    if verbose:
        print("\n########### PROCESSED WITH THE MEAN STELLAR ABUNDANCES ##############")

        single_simulation_report(planetabuN)


    # ============= MONTE-CARLO SIMULATION TO ESTIMATE PARAMETER DISTRIBUTION ===============#

    #### initialize output dictionaries
    mantle_keys = [
        "Na2O",
        "CaO",
        "MgO",
        "Al2O3",
        "SiO2",
        "FeO",
        "NiO",
        "SO3",
        "CO2",
        "C",
        "Metals",
        "ExtraO",
    ]
    core_keys = ["Fe", "Ni", "S", "Si"]

    mantle = dict.fromkeys(mantle_keys, [])
    mantle_Nval = dict.fromkeys(mantle_keys, [])
    core = dict.fromkeys(core_keys, [])
    core_Nval = dict.fromkeys(core_keys, [])

    fcoremass = np.zeros(0)

    # -------------- MC-loop -------------- #

    if verbose:
        print("############# MONTE-CARLO ###############")
        print("Performing "+str(Nsim)+" Monte-Carlo simulations using stellar abundance uncertainties ...")

    for n in tqdm(range(Nsim)):

        # draw random
        _planetabuN = planetabuN.copy()
        for elem in list(planetabuN.keys()):
            if not np.isnan(_planetabuN[elem]):

                # use asymmetric distribution
                _planetabuN[elem] = generate_asymmetric_random(
                    planetabuN[elem], planetabuNerrup[elem], planetabuNerrdn[elem]
                )

        # evaluate chemsysmodel on this random draw
        _mantle, _core, _fcoremass, Nvalc, _mantle_Nval, _core_Nval = chemsysmodel(
            _planetabuN, verbose=False
        )

        if Nvalc != 0 and Nvalc is not np.nan:
            for comp in mantle:
                mantle[comp] = np.append(mantle[comp], _mantle[comp])
                mantle_Nval[comp] = np.append(mantle_Nval[comp], _mantle_Nval[comp])
            for comp in core:
                core[comp] = np.append(core[comp], _core[comp])
                core_Nval[comp] = np.append(core_Nval[comp], _core_Nval[comp])
            fcoremass = np.append(fcoremass, _fcoremass)


    # ===================== GET OCCURANCE RATE ===================================#
    mantle_Nval_sum = np.nansum(list(mantle_Nval.values()), axis=1)
    core_Nval_sum = np.nansum(list(core_Nval.values()), axis=1)
    mantle_Nval_P = mantle_Nval_sum / sum(mantle_Nval["MgO"])
    core_Nval_P = core_Nval_sum / sum(mantle_Nval["MgO"])

    #breakpoint()
    # ===================== SCALE PARAMETERS WITH OCCURANCE RATE =====================#

    # ---- mantle ----#
    data = np.array(list(mantle.values()))
    data *= mantle_Nval_P[:, None]  # statistical weight accor. to occurance rate
    data /= np.nansum(data, axis=0)

    # get percentiles
    data = np.nanpercentile(data, q=[15.87, 50, 84.13], axis=1)

    data2 = data.copy()
    data2[0, :] = data[1, :] / np.nansum(data[1, :])  # *100.
    data2[1, :] = (data[2, :] - data[1, :]) / data[1, :] * data2[0, :]
    data2[2, :] = (data[1, :] - data[0, :]) / data[1, :] * data2[0, :]
    keys = np.array(list(mantle.keys()))
    data2 = np.vstack((np.array([keys]), data2))
    mantle_df = pandas.DataFrame(data2.T, columns=["param", "median", "errup", "errdn"])
   # mantle_df["N_val"] = mantle_Nval_P
    mantle_df = mantle_df.replace("nan", "0")

    mantle_df["median"] = pandas.to_numeric(mantle_df["median"]) * 100
    mantle_df["errup"] = pandas.to_numeric(mantle_df["errup"]) * 100
    mantle_df["errdn"] = pandas.to_numeric(mantle_df["errdn"]) * 100
    #mantle_df["N_val"] = pandas.to_numeric(mantle_df["N_val"])
    mantle_df.index = mantle_df["param"]

    #breakpoint()

    # ---- core ----#
    data = np.array(list(core.values()))
    data *= core_Nval_P[:, None]  # weight accor. to occurance rate
    data /= np.nansum(data, axis=0)

    # get percentiles
    data = np.nanpercentile(data, q=[15.87, 50, 84.13], axis=1)

    data2 = data.copy()
    data2[0, :] = data[1, :] / np.nansum(data[1, :])  # *100.
    data2[1, :] = (data[2, :] - data[1, :]) / data[1, :] * data2[0, :]
    data2[2, :] = (data[1, :] - data[0, :]) / data[1, :] * data2[0, :]
    keys = np.array(list(core.keys()))
    data2 = np.vstack((np.array([keys]), data2))
    core_df = pandas.DataFrame(data2.T, columns=["param", "median", "errup", "errdn"])
    #core_df["N_val"] = core_Nval_P
    core_df = core_df.replace("nan", "0")

    core_df["median"] = pandas.to_numeric(core_df["median"]) * 100
    core_df["errup"] = pandas.to_numeric(core_df["errup"]) * 100
    core_df["errdn"] = pandas.to_numeric(core_df["errdn"]) * 100
    #core_df["N_val"] = pandas.to_numeric(core_df["N_val"])
    core_df.index = core_df["param"]


    # ----- CMF -----#

    # get percentiles
    data = np.nanpercentile(fcoremass, q=[15.87, 50, 84.13])
    data2 = data.copy()
    data2[0] = data[1]
    data2[1] = data[2] - data[1]
    data2[2] = data[1] - data[0]
    data2 = np.append("CMF", data2)
    CMF_df = pandas.DataFrame(
        np.array([data2]), columns=["param", "median", "errup", "errdn"]
    )
    #CMF_df["N_val"] = 1.0

    CMF_df["median"] = pandas.to_numeric(CMF_df["median"]) * 100
    CMF_df["errup"] = pandas.to_numeric(CMF_df["errup"]) * 100
    CMF_df["errdn"] = pandas.to_numeric(CMF_df["errdn"]) * 100
    #CMF_df["N_val"] = pandas.to_numeric(CMF_df["N_val"])
    CMF_df.index = CMF_df["param"]

    #breakpoint()

    # ================= PRINT OUTPUT TO CONSOLE =======================#
    if verbose == True:
        print("\n\n############# FINAL RESULTS ###############\n")

        ExoInt_Output = ""

        ExoInt_Output += "------------- MANTLE (wt%) -------------\n"
        for param in mantle_df.param:
            ExoInt_Output += (
                param
                + ": "
                + str(round(mantle_df.loc[param]["median"], 2))
                + " + "
                + str(round(mantle_df.loc[param]["errup"], 2))
                + " - "
                + str(round(mantle_df.loc[param]["errdn"], 2))
                + "\n"
            )

        ExoInt_Output += "-------------- CORE (wt%) --------------\n"
        for param in core_df.param:
            ExoInt_Output += (
                param
                + ": "
                + str(round(core_df.loc[param]["median"], 2))
                + " + "
                + str(round(core_df.loc[param]["errup"], 2))
                + " - "
                + str(round(core_df.loc[param]["errdn"], 2))
                + "\n"
            )

        ExoInt_Output += "-------------- CORE (wt%) -------------- \n"
        ExoInt_Output += (
            "CMF"
            + ": "
            + str(round(CMF_df.loc["CMF"]["median"], 2))
            + " + "
            + str(round(CMF_df.loc["CMF"]["errup"], 2))
            + " - "
            + str(round(CMF_df.loc["CMF"]["errdn"], 2))
        )

        print(ExoInt_Output)

    # ================= WRITE OUTPUT TO FILE =======================#
    file = open(outdir + filename + "_exoE_results_formated.txt", "w")
    file.write(ExoInt_Output)
    file.close()

    results_df = pandas.concat([mantle_df, core_df, CMF_df])
    results_df.index = results_df["param"]
    del results_df["param"]
    results_df.to_csv(outdir + filename + "_exoE_results_all.csv")
    #breakpoint()
    ###output mantle, core, and cmf seperately, to be called for the subsequent mineral modellling     
    del mantle_df["param"], core_df["param"], CMF_df["param"]  ##remove the extra 'param' column for output
    mantle_df.to_csv(outdir + filename + '_exoE_mantlecomp_final.txt', sep='\t')
    core_df.to_csv(outdir + filename + '_exoE_corecomp_final.txt', sep='\t')
    CMF_dfc = CMF_df / 100. ##normalise to 1
    CMF_dfc.to_csv(outdir + filename + '_exoE_fcoremass_final.txt', sep='\t')
    #breakpoint()
    
    ##output estimated planetary abundance
    dl = {'elem':elemid, 'abu(Al=100)':np.array(planetabuN.values()), 'abu_errup':np.array(planetabuNerrup.values()), 'abu_errdn':np.array(planetabuNerrdn.values())}
    planetabu_df = pandas.DataFrame(dl)
    planetabu_df.to_csv(outdir + filename + '_exoE_bulkcomp_final.txt', sep='\t')
    # ===================== FINISH ===========================#

    print("\nDone.\n")

    return results_df
