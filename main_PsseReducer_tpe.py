# --------------------------------------------------------------------------------------------------
# Model reduction tool: Main function of two-port equivalent
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
    outpath = r"""TpeOut\\"""

    # TPE reduction
    TPE_list = r"""TPE_list.xlsx"""
    rootbus1_list, bus_red_list, d1bus_list, rootbus2_list, retdbus_list, N_tpe = read_bus_TPE(TPE_list)

    for instance_i in range(N_tpe):
        bus_root1 = rootbus1_list[instance_i]
        bus_root2 = rootbus2_list[instance_i]
        bus_ret = retdbus_list[instance_i]
        bus_red = bus_red_list[instance_i]
        bus_d1 = d1bus_list[instance_i]

        # read power flow model before reduction step i
        psspy.psseinit(50000)
        if instance_i == 0:
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


        # prepare data for two-port equivalent
        P1, Q1, Vm1, Va1, PrateA1, PrateB1, P3, Q3, Vm3, Va3, PrateA3, PrateB3, PL, QL, PG, QG, MW_ll, MW_ul, Mvar_ll,\
        Mvar_ul, MVA_base, PS, QS, load_bus, gen_bus, shunt_bus = read_subsys_TPE(pfd, bus_root1, bus_red, bus_d1, bus_root2, bus_ret)



        # calculate TPE equivalent
        Vm2, Va2, r1, x1, r2, x2 = CalcTwoPortEqui(pfd, P1, Q1, Vm1, Va1, P3, Q3, Vm3, Va3, PL, QL, PG, QG, PS, QS)



        # Implement two-port equivalent in PSSE
        DoTpeInPsse(psspy, pfd, bus_root1, bus_root2, bus_ret, bus_red, bus_d1, PL, QL, PG, QG, PS, QS, PrateA1, PrateB1, PrateA3,
                PrateB3, MW_ll, MW_ul, Mvar_ll, Mvar_ul, MVA_base, Vm2, Va2, r1, x1, r2, x2)



        # save new power flow data
        psspy.rawd_2(0, 1, [0, 0, 1, 0, 0, 0, 0], 0, outpath + "tpe_step_" + str(instance_i + 1))



    # calc a summary for TPE reduction
    pfd.getdata(psspy)
    n_gen_af, n_load_af, n_bus_af, n_line_af, n_xfmr_af, n_shunt_af = CountEle(pfd)



    #
    print("\nTPE reduction summary:")
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








