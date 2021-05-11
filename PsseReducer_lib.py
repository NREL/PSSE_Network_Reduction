# --------------------------------------------------------------------------------------------------
# Model reduction tool functions
# Bin Wang
# 5/11/2021
#
# Reference:
# Wang, Bin, Andy Hoke, and Jin Tan. 2021. Power System Network Reduction for Power Hardware-in-the-Loop Simulation: Preprint.
# Golden, CO: National Renewable Energy Laboratory. NREL/CP-5D00-78372. https://www.nrel.gov/docs/fy21osti/78372.pdf
# --------------------------------------------------------------------------------------------------

import xlrd
import numpy
import math


# --------------------------------------------------------------------------------------------------
# power flow data class
class PFData():
    def __init__(self):
        # system data
        self.basemva = []

        # bus data
        self.bus_num = []
        self.bus_type = []
        self.bus_Vm = []
        self.bus_Va = []
        self.bus_kV = []
        self.bus_basekV = []
        self.bus_name = []
        self.bus_area = []
        self.bus_zone = []
        self.bus_owner = []

        # load data
        self.load_id = []
        self.load_bus = []
        self.load_Z = []
        self.load_I = []
        self.load_P = []
        self.load_MW = []
        self.load_Mvar = []
        self.load_MW_total = []
        self.load_Mvar_total = []

        # generator data
        self.gen_id = []
        self.gen_bus = []
        self.gen_S = []
        self.gen_mod = []
        self.gen_MW = []
        self.gen_Mvar = []
        self.gen_MW_total = []
        self.gen_Mvar_total = []
        self.gen_MW_ll = []
        self.gen_MW_ul = []
        self.gen_Mvar_ll = []
        self.gen_Mvar_ul = []
        self.gen_MVA_base = []

        # branch data
        self.brc_from = []
        self.brc_to = []
        self.brc_id = []
        self.brc_rateA = []
        self.brc_rateB = []
        self.brc_S = []
        self.brc_P = []
        self.brc_Q = []
        self.brc_line_id = []
        self.brc_xfmr_id = []

        # shunt data
        self.shunt_id = []
        self.shunt_bus = []
        self.shunt_S_NOM = []
        self.shunt_MW_NOM = []
        self.shunt_Mvar_NOM = []


    def getdata(self, psspy):
        # system data
        self.basemva = psspy.sysmva()

        # bus data
        self.bus_num = psspy.abusint(-1, 2, 'NUMBER')[1][0]
        self.bus_type = psspy.abusint(-1, 2, 'TYPE')[1][0]
        self.bus_Vm = psspy.abusreal(-1, 2, 'PU')[1][0]
        self.bus_Va = psspy.abusreal(-1, 2, 'ANGLE')[1][0]
        self.bus_kV = psspy.abusreal(-1, 2, 'KV')[1][0]
        self.bus_basekV = psspy.abusreal(-1, 2, 'BASE')[1][0]
        self.bus_name = psspy.abuschar(-1, 1, 'NAME')[1][0]
        self.bus_area = psspy.abusint(-1, 1, 'AREA')[1][0]
        self.bus_zone = psspy.abusint(-1, 1, 'ZONE')[1][0]
        self.bus_owner = psspy.abusint(-1, 1, 'OWNER')[1][0]

        # load data
        self.load_id = psspy.aloadchar(-1, 1, 'ID')[1][0]
        self.load_bus = psspy.aloadint(-1, 1, 'NUMBER')[1][0]
        self.load_Z = numpy.asarray(psspy.aloadcplx(-1, 1, 'YLACT')[1][0])
        self.load_I = numpy.asarray(psspy.aloadcplx(-1, 1, 'ILACT')[1][0])
        self.load_P = numpy.asarray(psspy.aloadcplx(-1, 1, 'MVAACT')[1][0])
        self.load_MW = self.load_Z.real + self.load_I.real + self.load_P.real
        self.load_Mvar = self.load_Z.imag + self.load_I.imag + self.load_P.imag
        self.load_MW_total = sum(self.load_MW)
        self.load_Mvar_total = sum(self.load_Mvar)

        # generator data
        self.gen_id = psspy.amachchar(-1, 1, 'ID')[1][0]
        self.gen_bus = psspy.amachint(-1, 1, 'NUMBER')[1][0]
        self.gen_S = numpy.asarray(psspy.amachcplx(-1, 1, 'PQGEN')[1][0])
        self.gen_mod = numpy.asarray(psspy.amachint(-1, 1, 'WMOD')[1][0])
        self.gen_MW = self.gen_S.real
        self.gen_Mvar = self.gen_S.imag
        self.gen_MW_total = sum(self.gen_MW)
        self.gen_Mvar_total = sum(self.gen_Mvar)
        self.gen_MW_ll = psspy.amachreal(-1, 1, 'PMIN')[1][0]
        self.gen_MW_ul = psspy.amachreal(-1, 1, 'PMAX')[1][0]
        self.gen_Mvar_ll = psspy.amachreal(-1, 1, 'QMIN')[1][0]
        self.gen_Mvar_ul = psspy.amachreal(-1, 1, 'QMAX')[1][0]
        self.gen_MVA_base = psspy.amachreal(-1, 1, 'MBASE')[1][0]

        # branch data
        ierr, iarray = psspy.abrnint(-1, 0, 0, 3, 2, ['FROMNUMBER', 'TONUMBER'])
        self.brc_from = iarray[0][:]
        self.brc_to = iarray[1][:]
        self.brc_id = psspy.abrnchar(-1, 0, 0, 3, 2, ['ID'])[1][0]
        self.brc_rateA = psspy.abrnreal(-1, 0, 0, 3, 2, ['RATEA'])[1][0]
        self.brc_rateB = psspy.abrnreal(-1, 0, 0, 3, 2, ['RATEB'])[1][0]
        self.brc_S = numpy.asarray(psspy.abrncplx(-1, 1, 1, 3, 2, ['PQ'])[1][0])
        self.brc_P = self.brc_S.real
        self.brc_Q = self.brc_S.imag
        self.brc_line_id = psspy.abrnchar(-1, 0, 0, 1, 1, ['ID'])[1][0]
        self.brc_xfmr_id = psspy.abrnchar(-1, 0, 0, 5, 1, ['ID'])[1][0]

        # shunt data
        self.shunt_id = psspy.afxshuntchar(-1, 1, 'ID')[1][0]
        self.shunt_bus = psspy.afxshuntint(-1, 1, 'NUMBER')[1][0]
        self.shunt_S_NOM = numpy.asarray(psspy.afxshuntcplx(-1, 1, 'SHUNTNOM')[1][0])
        self.shunt_MW_NOM = self.shunt_S_NOM.real
        self.shunt_Mvar_NOM = self.shunt_S_NOM.imag




def read_bus(SPE_list):
    list_file = xlrd.open_workbook(SPE_list)
    sheet = list_file.sheet_by_index(0)
    N_list = sheet.ncols

    rootbus_list = []
    d1bus_list = []
    redbus_list = []

    for coli in range(N_list):
        N_len_i = int(sheet.cell_value(0, coli))
        N_d1bus_i = int(sheet.cell_value(1, coli))

        rootbus_i = int(sheet.cell_value(2, coli))

        d1bus_i = []
        for ii in range(N_d1bus_i):
            tempbus = sheet.cell_value(ii + 3, coli)
            if tempbus:
                d1bus_i.append(int(tempbus))

        redbus_i = []
        for ii in range(N_len_i - 1):
            redbus_i.append(int(sheet.cell_value(ii + 3, coli)))

        rootbus_list.append(rootbus_i)
        d1bus_list.append(d1bus_i)
        redbus_list.append(redbus_i)



    return rootbus_list, redbus_list, d1bus_list, N_list
    


def read_subsys(pfd, bus_rut, bus_red, bus_d1):
    # combine all loads on buses in reduced area
    PL = 0  # MW
    QL = 0  # Mvar
    load_bus = []
    for i in range(len(pfd.load_bus)):      # consider multiple loads on the same bus
        busi = pfd.load_bus[i]
        if busi in bus_red:
            PL = PL + pfd.load_MW[i]
            QL = QL + pfd.load_Mvar[i]
            if busi not in load_bus:
                load_bus.append(busi)

    # combine all shunts on buses in reduced area
    PS = 0.0  # MW
    QS = 0.0  # Mvar
    shunt_bus = []
    for i in range(len(pfd.shunt_bus)):  # consider multiple shunts on the same bus
        busi = pfd.shunt_bus[i]
        busi_idx = pfd.bus_num.index(busi)
        if busi in bus_red:
            PS = PS + pfd.shunt_MW_NOM[i]
            QS = QS + pfd.shunt_Mvar_NOM[i]
            if busi not in shunt_bus:
                shunt_bus.append(busi)


    # combine all generations on buses in reduced area
    PG = 0.0  # MW
    QG = 0.0  # MW

    MW_ll = 0.0
    MW_ul = 0.0
    Mvar_ll = 0.0
    Mvar_ul = 0.0
    MVA_base = 0.0

    gen_bus = []
    for busi in bus_red:
        last_found = -1
        while busi in pfd.gen_bus[last_found + 1:]:  # consider multiple gens on the same bus
            last_found = pfd.gen_bus.index(busi, last_found + 1)
            if last_found == -1:
                break
            PG = PG + pfd.gen_MW[last_found]
            QG = QG + pfd.gen_Mvar[last_found]

            MW_ll = MW_ll + pfd.gen_MW_ll[last_found]
            MW_ul = MW_ul + pfd.gen_MW_ul[last_found]
            Mvar_ll = Mvar_ll + pfd.gen_Mvar_ll[last_found]
            Mvar_ul = Mvar_ul + pfd.gen_Mvar_ul[last_found]
            MVA_base = MVA_base + pfd.gen_MVA_base[last_found]

            if busi not in gen_bus:
                gen_bus.append(busi)



    # get the average voltage magnitude
    n = 0
    vn = 0
    for busi in bus_red:
        last_found = -1
        while busi in pfd.load_bus[last_found + 1:]:
            last_found = pfd.load_bus.index(busi, last_found + 1)
            if last_found == -1:
                break
            n = n + 1
            bus_idx = pfd.bus_num.index(busi)
            vn = vn + pfd.bus_Vm[bus_idx]

        last_found = -1
        while busi in pfd.gen_bus[last_found + 1:]:
            last_found = pfd.gen_bus.index(busi, last_found + 1)
            if last_found == -1:
                break
            n = n + 1
            bus_idx = pfd.bus_num.index(busi)
            vn = vn + pfd.bus_Vm[bus_idx]
    Ve = vn / n



    # get Pin, Qin, Vm and Va at root bus
    bus_idx = pfd.bus_num.index(bus_rut)
    Vm = pfd.bus_Vm[bus_idx]
    Va = pfd.bus_Va[bus_idx]

    Pin = 0
    Qin = 0
    last_found = -1
    PrateA = 0
    PrateB = 0
    while bus_rut in pfd.brc_from[last_found + 1:]:
        last_found = pfd.brc_from.index(bus_rut, last_found + 1)
        if pfd.brc_to[last_found] in bus_d1:
            Pin = Pin + pfd.brc_P[last_found]
            Qin = Qin + pfd.brc_Q[last_found]
            PrateA = PrateA + pfd.brc_rateA[last_found]
            PrateB = PrateB + pfd.brc_rateB[last_found]

    return Pin, Qin, PL, QL, PG, QG, PS, QS, Vm, Va, Ve, PrateA, PrateB, MW_ll, MW_ul, Mvar_ll, Mvar_ul, MVA_base, load_bus, gen_bus, shunt_bus



def CalcSinglePortEqui(pfd, Pin, Qin, PL, QL, PG, QG, PS, QS, Vm, Va, Ve):
    Sin = numpy.complex(Pin, Qin) / pfd.basemva
    V1 = numpy.complex(Vm * math.cos(Va), Vm * math.sin(Va))
    alpha = Sin / V1 * Ve
    Se = numpy.complex(PL - PG - PS, QL - QG - QS) / pfd.basemva
    a = alpha.real
    b = alpha.imag
    c = Se.real
    d = Se.imag

    k = Ve / Vm * numpy.abs(Sin) / numpy.abs(Se)
    Vae = math.asin(k * d / math.sqrt(a * a + b * b)) - math.atan(b / a)

    V2 = numpy.complex(Ve * math.cos(Vae), Ve * math.sin(Vae))
    Z = (V1 - V2/k)/numpy.conj(Sin/V1)
    r = Z.real
    x = Z.imag
    return k, Vae, r, x







def DoSpeInPsse(psspy, bus_root, bus_red, bus_d1, pfd, PL, QL, PG, QG, PS, QS, Ve, PrateA, PrateB, MW_ll, MW_ul, Mvar_ll, Mvar_ul, MVA_base, k, r, x):
    Ve_kV = 20.0  # set kV at equivalent bus

    # implement equivalent
    # remove all lines/transformers in reduced area
    bus_allred = bus_red
    bus_allred.append(bus_root)
    for brc_from, brc_to, brc_id in zip(pfd.brc_from, pfd.brc_to, pfd.brc_id):
        if brc_from in bus_allred:
            if brc_to in bus_allred:
                psspy.purgbrn(brc_from, brc_to, brc_id)

    # remove all gen in reduced area
    bus_allred.remove(bus_root)
    for gen_bus, gen_id in zip(pfd.gen_bus, pfd.gen_id):
        if gen_bus in bus_allred:
            psspy.purgmac(gen_bus, gen_id)

    # remove all loads in reduced area
    for load_bus, load_id in zip(pfd.load_bus, pfd.load_id):
        if load_bus in bus_allred:
            psspy.purgload(load_bus, load_id)

    # remove all shunts in reduced area
    for shunt_bus, shunt_id in zip(pfd.shunt_bus, pfd.shunt_id):
        if shunt_bus in bus_allred:
            psspy.purgshunt(shunt_bus, shunt_id)


    # remove all buses except for bus_root and bus_d1[0] in reduced area
    bus_removed = bus_red
    bus_removed.remove(bus_d1[0])
    for busi in bus_removed:
        psspy.bsysinit(1)
        psspy.bsyso(1, int(busi))
        psspy.extr(1, 0, [0, 0])

    # change voltage level and bus type to PV if PG>0
    if PG > 0:
        # if reduced area contains original slack bus, then change to Slack bus, otherwise, to a PV bus
        SL = pfd.bus_type.index(3)
        if pfd.bus_num[SL] in bus_removed:
            psspy.bus_chng_3(bus_d1[0], [3, psspy._i, psspy._i, psspy._i],
                             [Ve_kV, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f], psspy._s)
        else:
            psspy.bus_chng_3(bus_d1[0], [2, psspy._i, psspy._i, psspy._i],
                             [Ve_kV, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f], psspy._s)
    else:
        psspy.bus_chng_3(bus_d1[0], [1, psspy._i, psspy._i, psspy._i],
                         [Ve_kV, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f], psspy._s)


    # add transformer
    psspy.two_winding_data_4(bus_root, bus_d1[0], r"""1""",
                             [psspy._i, bus_d1[0], psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i,
                              bus_d1[0],
                              psspy._i, psspy._i, 0, psspy._i, psspy._i, psspy._i],
                             [r, x, psspy._f, k, psspy._f, psspy._f, psspy._f, psspy._f, PrateA, PrateB, psspy._f,
                              psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
                              psspy._f,
                              psspy._f, psspy._f,
                              psspy._f, psspy._f], ["", r""""""])



    # add gen
    if PG > 0:
        psspy.plant_data(bus_d1[0], psspy._i, [Ve, psspy._f])
        psspy.machine_data_2(bus_d1[0], r"""eq""", [psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, 1],
                             [PG, psspy._f, Mvar_ul, Mvar_ll, MW_ul, MW_ll, MVA_base, 0.0,
                              0.25, psspy._f, psspy._f, psspy._f,
                              psspy._f, psspy._f, psspy._f, psspy._f, psspy._f])



    # add load
    psspy.load_data_4(bus_d1[0], r"""eq""", [psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i],
                      [PL, QL, psspy._f, psspy._f, psspy._f, psspy._f])

    # add shunt
    if math.fabs(QS) > 0:
        psspy.shunt_data(bus_d1[0], r"""eq""", psspy._i, [PS, QS])

    # solve power flow
    psspy.fnsl([1, 0, 0, 1, 1, 0, 0, 0])
    psspy.fnsl([1, 0, 0, 1, 1, 0, 0, 0])



def CountEle(pfd):
    n_gen = len(pfd.gen_id)
    n_load = len(pfd.load_bus)
    n_bus = len(pfd.bus_num)
    n_line = len(pfd.brc_line_id)
    n_xfmr = len(pfd.brc_xfmr_id)
    n_shunt = len(pfd.shunt_id)
    return n_gen, n_load, n_bus, n_line, n_xfmr, n_shunt








def read_bus_TPE(TPE_list):
    list_file = xlrd.open_workbook(TPE_list)
    sheet = list_file.sheet_by_index(0)
    N_list = sheet.ncols

    rootbus1_list = []
    redbus_list = []
    d1bus_list = []
    rootbus2_list = []
    retdbus_list = []

    for coli in range(N_list):
        N_len_i = int(sheet.cell_value(0, coli))
        N_d1bus_i = int(sheet.cell_value(1, coli))

        rootbus1_i = int(sheet.cell_value(2, coli))
        rootbus2_i = int(sheet.cell_value(N_len_i + 1, coli))
        retdbus_i = int(sheet.cell_value(N_len_i, coli))

        d1bus_i = []
        for ii in range(N_d1bus_i):
            tempbus = sheet.cell_value(ii + 3, coli)
            if tempbus:
                d1bus_i.append(int(tempbus))

        redbus_i = []
        for ii in range(N_len_i - 2):
            redbus_i.append(int(sheet.cell_value(ii + 3, coli)))

        rootbus1_list.append(rootbus1_i)
        rootbus2_list.append(rootbus2_i)
        retdbus_list.append(retdbus_i)
        d1bus_list.append(d1bus_i)
        redbus_list.append(redbus_i)


    return rootbus1_list, redbus_list, d1bus_list, rootbus2_list, retdbus_list, N_list




def read_subsys_TPE(pfd, bus_root1, bus_red, bus_d1, bus_root2, bus_ret):
    # combine all loads on buses in reduced area
    PL = 0.0  # MW
    QL = 0.0  # Mvar
    load_bus = []
    for i in range(len(pfd.load_bus)):      # consider multiple loads on the same bus
        busi = pfd.load_bus[i]
        if busi in bus_red:
            PL = PL + pfd.load_MW[i]
            QL = QL + pfd.load_Mvar[i]
            if busi not in load_bus:
                load_bus.append(busi)

    # combine all shunts on buses in reduced area
    PS = 0.0  # MW
    QS = 0.0  # Mvar
    shunt_bus = []
    for i in range(len(pfd.shunt_bus)):  # consider multiple shunts on the same bus
        busi = pfd.shunt_bus[i]
        busi_idx = pfd.bus_num.index(busi)
        if busi in bus_red:
            PS = PS + pfd.shunt_MW_NOM[i]
            QS = QS + pfd.shunt_Mvar_NOM[i]
            if busi not in shunt_bus:
                shunt_bus.append(busi)


    # get injection from bus_root2
    P3 = 0.0
    Q3 = 0.0
    PrateA3 = 0.0
    PrateB3 = 0.0
    last_found = -1
    while bus_ret in pfd.brc_to[last_found + 1:]:
        last_found = pfd.brc_to.index(bus_ret, last_found + 1)
        if pfd.brc_from[last_found] == bus_root2:
            P3 = P3 + pfd.brc_P[last_found]
            Q3 = Q3 + pfd.brc_Q[last_found]
            PrateA3 = PrateA3 + pfd.brc_rateA[last_found]
            PrateB3 = PrateB3 + pfd.brc_rateB[last_found]

    bus_idx = pfd.bus_num.index(bus_root2)
    Vm3 = pfd.bus_Vm[bus_idx]
    Va3 = pfd.bus_Va[bus_idx]


    # combine all generations on buses in reduced area
    PG = 0.0  # MW
    QG = 0.0  # MW

    MW_ll = 0.0
    MW_ul = 0.0
    Mvar_ll = 0.0
    Mvar_ul = 0.0
    MVA_base = 0.0


    gen_bus = []
    for busi in bus_red:
        last_found = -1
        while busi in pfd.gen_bus[last_found + 1:]:  # consider multiple gens on the same bus
            last_found = pfd.gen_bus.index(busi, last_found + 1)
            if last_found == -1:
                break
            PG = PG + pfd.gen_MW[last_found]
            QG = QG + pfd.gen_Mvar[last_found]

            MW_ll = MW_ll + pfd.gen_MW_ll[last_found]
            MW_ul = MW_ul + pfd.gen_MW_ul[last_found]
            Mvar_ll = Mvar_ll + pfd.gen_Mvar_ll[last_found]
            Mvar_ul = Mvar_ul + pfd.gen_Mvar_ul[last_found]
            MVA_base = MVA_base + pfd.gen_MVA_base[last_found]

            if busi not in gen_bus:
                gen_bus.append(busi)


    # get Pin, Qin, Vm and Va at root bus1
    bus_idx = pfd.bus_num.index(bus_root1)
    Vm1 = pfd.bus_Vm[bus_idx]
    Va1 = pfd.bus_Va[bus_idx]

    P1 = 0.0
    Q1 = 0.0
    last_found = -1
    PrateA1 = 0.0
    PrateB1 = 0.0
    while bus_root1 in pfd.brc_from[last_found + 1:]:
        last_found = pfd.brc_from.index(bus_root1, last_found + 1)
        if pfd.brc_to[last_found] in bus_d1:
            P1 = P1 + pfd.brc_P[last_found]
            Q1 = Q1 + pfd.brc_Q[last_found]
            PrateA1 = PrateA1 + pfd.brc_rateA[last_found]
            PrateB1 = PrateB1 + pfd.brc_rateB[last_found]


    return P1, Q1, Vm1, Va1, PrateA1, PrateB1, P3, Q3, Vm3, Va3, PrateA3, PrateB3, PL, QL, PG, QG, MW_ll, MW_ul, Mvar_ll, Mvar_ul, MVA_base, PS, QS, load_bus, gen_bus, shunt_bus




def CalcTwoPortEqui(pfd, P1, Q1, Vm1, Va1, P3, Q3, Vm3, Va3, PL, QL, PG, QG, PS, QS):
    Sin1 = numpy.complex(P1, Q1) / pfd.basemva
    Sin3 = numpy.complex(P3, Q3) / pfd.basemva
    Se = numpy.complex(PG + PS - PL, QG + QS - QL) / pfd.basemva
    V1 = numpy.complex(Vm1 * math.cos(Va1), Vm1 * math.sin(Va1))
    V3 = numpy.complex(Vm3 * math.cos(Va3), Vm3 * math.sin(Va3))

    if numpy.abs(Se) < 1e-5:
        V2 = (V1+V3)/2
        Vm2 = numpy.abs(V2)
        Va2 = numpy.angle(V2)
    else:
        V2 = - Se / (Sin1 / V1 + Sin3 / V3) # eq.12 in the paper lacks a minus sign
        Vm2 = numpy.abs(V2)
        Va2 = numpy.angle(V2)
    Z1 = (V1 - V2) * numpy.conj(V1) / numpy.conj(Sin1)
    Z2 = (V3 - V2) * numpy.conj(V3) / numpy.conj(Sin3)

    r1 = Z1.real
    x1 = Z1.imag
    r2 = Z2.real
    x2 = Z2.imag
    return Vm2, Va2, r1, x1, r2, x2




def DoTpeInPsse(psspy, pfd, bus_root1, bus_root2, bus_ret, bus_red, bus_d1, PL, QL, PG, QG, PS, QS, PrateA1, PrateB1,
                PrateA3, PrateB3, MW_ll, MW_ul, Mvar_ll, Mvar_ul, MVA_base, Vm2, Va2, r1, x1, r2, x2):
    # implement equivalent
    # remove all lines/transformers in reduced area
    bus_allred = bus_red
    bus_allred.append(bus_root1)
    bus_allred.append(bus_root2)
    for brc_from, brc_to, brc_id in zip(pfd.brc_from, pfd.brc_to, pfd.brc_id):
        if brc_from in bus_allred:
            if brc_to in bus_allred:
                psspy.purgbrn(brc_from, brc_to, brc_id)

    # remove all gen in reduced area
    bus_allred.remove(bus_root1)
    bus_allred.remove(bus_root2)
    for gen_bus, gen_id in zip(pfd.gen_bus, pfd.gen_id):
        if gen_bus in bus_allred:
            psspy.purgmac(gen_bus, gen_id)

    # remove all loads in reduced area
    for load_bus, load_id in zip(pfd.load_bus, pfd.load_id):
        if load_bus in bus_allred:
            psspy.purgload(load_bus, load_id)

    # remove all shunts in reduced area
    for shunt_bus, shunt_id in zip(pfd.shunt_bus, pfd.shunt_id):
        if shunt_bus in bus_allred:
            psspy.purgshunt(shunt_bus, shunt_id)


    # remove all buses except for bus_root and bus_ret in reduced area
    bus_removed = bus_red
    bus_removed.remove(bus_ret)
    for busi in bus_removed:
        psspy.bsysinit(1)
        psspy.bsyso(1, int(busi))
        psspy.extr(1, 0, [0, 0])

    # change bus type, bus voltage mag and angle at retained bus
    idx = pfd.bus_num.index(bus_root1)
    kV_root = pfd.bus_basekV[idx]
    if PG > 0:
        # if reduced area contains original slack bus, then change to Slack bus, otherwise, to a PV bus
        SL = pfd.bus_type.index(3)
        if pfd.bus_num[SL] in bus_removed:
            psspy.bus_chng_3(bus_ret, [3, psspy._i, psspy._i, psspy._i],
                             [kV_root, Vm2, Va2*180/3.1415926, psspy._f, psspy._f, psspy._f, psspy._f], psspy._s)
        else:
            psspy.bus_chng_3(bus_ret, [2, psspy._i, psspy._i, psspy._i],
                             [kV_root, Vm2, Va2*180/3.1415926, psspy._f, psspy._f, psspy._f, psspy._f], psspy._s)
    else:
        psspy.bus_chng_3(bus_ret, [1, psspy._i, psspy._i, psspy._i],
                         [kV_root, Vm2, Va2*180/3.1415926, psspy._f, psspy._f, psspy._f, psspy._f], psspy._s)


    # add line
    psspy.branch_data(bus_root1, bus_ret, r"""eq""", [psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i],
                      [r1, x1, psspy._f, PrateA1, PrateB1, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
                       psspy._f, psspy._f, psspy._f, psspy._f])
    psspy.branch_data(bus_root2, bus_ret, r"""eq""", [psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i],
                      [r2, x2, psspy._f, PrateA3, PrateB3, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f, psspy._f,
                       psspy._f, psspy._f, psspy._f, psspy._f])



    # add gen
    if PG > 0:
        psspy.plant_data(bus_ret, psspy._i, [Vm2, psspy._f])
        psspy.machine_data_2(bus_ret, r"""eq""", [psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, 1],
                             [PG, psspy._f, Mvar_ul, Mvar_ll, MW_ul, MW_ll, MVA_base, 0.0,
                              0.25, psspy._f, psspy._f, psspy._f,
                              psspy._f, psspy._f, psspy._f, psspy._f, psspy._f])



    # add load
    if math.fabs(PL) > 0:
        psspy.load_data_4(bus_ret, r"""eq""", [psspy._i, psspy._i, psspy._i, psspy._i, psspy._i, psspy._i],
                          [PL, QL, psspy._f, psspy._f, psspy._f, psspy._f])


    # add shunt
    if math.fabs(QS) > 0:
        psspy.shunt_data(bus_ret, r"""eq""", psspy._i, [PS, QS])


    # psspy.rawd_2(0, 1, [0, 0, 1, 0, 0, 0, 0], 0, "TpeOut\\tpe_step_temp")

    # solve power flow
    psspy.fnsl([1, 0, 0, 1, 1, 0, 0, 0])
    psspy.fnsl([1, 0, 0, 1, 1, 0, 0, 0])
