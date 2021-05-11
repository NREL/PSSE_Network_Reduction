# --------------------------------------------------------------------------------------------------
# Model reduction tool: Main function of single-port equivalent
# Bin Wang
# 5/11/2021
#
# Reference:
# Wang, Bin, Andy Hoke, and Jin Tan. 2021. Power System Network Reduction for Power Hardware-in-the-Loop Simulation: Preprint.
# Golden, CO: National Renewable Energy Laboratory. NREL/CP-5D00-78372. https://www.nrel.gov/docs/fy21osti/78372.pdf
# --------------------------------------------------------------------------------------------------


import os, sys

sys_path_PSSE = r'C:\Program Files (x86)\PTI\PSSE33\PSSBIN'
sys.path.append(sys_path_PSSE)
os_path_PSSE = r'C:\Program Files (x86)\PTI\PSSE33\PSSBIN'
os.environ['PATH'] += ';' + os_path_PSSE
os.environ['PATH'] += ';' + sys_path_PSSE

import psspy
from PsseReducer_lib import *




def main():
    # set output path
    outpath = r"""SpeOut\\"""

    # SPE reduction
    SPE_list = r"""SPE_list.xlsx"""
    rootbus_list, redbus_list, d1bus_list, N_spe = read_bus(SPE_list)

    for instance_i in range(N_spe):
        bus_root = rootbus_list[instance_i]
        bus_red = redbus_list[instance_i]
        bus_d1 = d1bus_list[instance_i]

        # read power flow model before reduction step i
        psspy.psseinit(50000)
        if instance_i == 0:
            # psspy.readrawversion(0, r"""30""", r"""FullModel\wecc179.raw""")  # reading .raw file
            # psspy.case(r"""FullModel\Maui2022dm_v4_wHydro_step0""")  # reading .sav file
            psspy.read(0, r"""FullModel\wecc179_v33.raw""")
        else:
            pass
        psspy.fnsl([1, 0, 0, 1, 1, 0, 0, 0])

        # get power flow data
        pfd = PFData()
        pfd.getdata(psspy)

        # count elements in full model
        if instance_i == 0:
            n_gen_bf, n_load_bf, n_bus_bf, n_line_bf, n_xfmr_bf, n_shunt_bf = CountEle(pfd)


        # prepare data for single-port equivalent
        Pin, Qin, PL, QL, PG, QG, PS, QS, Vm, Va, Ve, PrateA, PrateB, MW_ll, MW_ul, Mvar_ll, Mvar_ul, MVA_base, \
        red_load_bus, red_gen_bus, red_shunt_bus = read_subsys(pfd, bus_root, bus_red, bus_d1)



        # calculate SPE equivalent
        k, Vae, r, x = CalcSinglePortEqui(pfd, Pin, Qin, PL, QL, PG, QG, PS, QS, Vm, Va, Ve)



        # Implement single-port equivalent in PSSE
        DoSpeInPsse(psspy, bus_root, bus_red, bus_d1, pfd, PL, QL, PG, QG, PS, QS, Ve, PrateA, PrateB, MW_ll, MW_ul,
                       Mvar_ll, Mvar_ul, MVA_base, k, r, x)



        # save new power flow data
        psspy.rawd_2(0, 1, [0, 0, 1, 0, 0, 0, 0], 0, outpath + "spe_step_" + str(instance_i + 1))



    # calc a summary for SPE reduction
    pfd.getdata(psspy)
    n_gen_af, n_load_af, n_bus_af, n_line_af, n_xfmr_af, n_shunt_af = CountEle(pfd)



    #
    print("\nSPE reduction summary:")
    print("(# of elements, Before, After, Reduction %)")
    print("-----------------------------------------------")
    print("       Buses: ", n_bus_bf, n_bus_af, str(float((n_bus_bf - n_bus_af))/n_bus_bf*100)[0:5] + "%")
    print(" Generations: ", n_gen_bf, n_gen_af, str(float((n_gen_bf - n_gen_af))/n_gen_bf*100)[0:5] + "%")
    print("       Loads: ", n_load_bf, n_load_af, str(float((n_load_bf - n_load_af))/n_load_bf*100)[0:5] + "%")
    print("       Lines: ", n_line_bf, n_line_af, str(float((n_line_bf - n_line_af))/n_line_bf*100)[0:5] + "%")
    print("Transformers: ", n_xfmr_bf + 1, n_xfmr_af + 1, str(float((n_xfmr_bf - n_xfmr_af))/(n_xfmr_bf+1)*100)[0:5] + "%")
    print("      Shunts: ", n_shunt_bf, n_shunt_af,
          str(float((n_shunt_bf - n_shunt_af)) / (n_shunt_bf + 1) * 100)[0:5] + "%")




# main function
main()








