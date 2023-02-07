# Force Field to calculate sugar
# Bond distance, angles, Moorse potential(for HB), partial charges
# All the 4 calculated energies can be added together and compared with another structure of interest
# This file can already be used on multiple structures
# Most updated: 24/10/2021
from csv_tools import read_csv_file
from os_tools import file_list
from sugar_tools import Sugar
from os.path import basename
from os_tools import single_cmd
import numpy
import math
import pandas as pd
from gaussian_tools import out2zpe

mm_ls_type = read_csv_file("bonding_data.csv")['000_TEXT']
mm_ls_min = read_csv_file("bonding_data.csv")['001_TEXT']
mm_ls_max = read_csv_file("bonding_data.csv")['002_TEXT']

bond_para = read_csv_file("force_constants_bond.csv")
typeb = list(bond_para['000_TEXT'])
xeqb = list(bond_para['001_TEXT'])
xconsb = []
xconsb.append(list(bond_para['002_TEXT'])) #x4
xconsb.append(list(bond_para['003_TEXT'])) #x3
xconsb.append(list(bond_para['004_TEXT'])) #x2
xconsb.append(list(bond_para['005_TEXT'])) #x1
xconsb.append(list(bond_para['006_TEXT'])) #constant


angle_para = read_csv_file("force_constants_angle.csv")
typea = list(angle_para['000_TEXT'])
xeqa = list(angle_para['001_TEXT'])
xconsa = []
xconsa.append(list(angle_para['002_TEXT'])) #x4
xconsa.append(list(angle_para['003_TEXT'])) #x3
xconsa.append(list(angle_para['004_TEXT'])) #x2
xconsa.append(list(angle_para['005_TEXT'])) #x1
xconsa.append(list(angle_para['006_TEXT'])) #constant


type_ls = []
min_ls = []
max_ls = []
for h in range(1, len(mm_ls_type)):
    type_ls.append(mm_ls_type[h])
    min_ls.append(mm_ls_min[h])
    max_ls.append(mm_ls_max[h])

def dis_func(bb,coor):  # Collect all bonds and showed in a list with [distance, bond_type, atom_idx]
    #coor = Sugar(bb).next()
    dis = []
    dis_t = []
    dis_idx = []
    for i in range(len(coor)):
        ai = numpy.array(coor[i][1:])
        ait = coor[i][0]
        for j in range(i+1, len(coor)):
            aj = numpy.array(coor[j][1:])
            cur_dis = numpy.linalg.norm(ai - aj)
            ajt = coor[j][0]
            bond_t = ait + '-' + ajt
            dis.append(cur_dis)
            dis_t.append(bond_t)
            ip = i + 1
            jp = j + 1
            dis_idx.append(str(ip)+'-'+str(jp))

    useful_dis = []
    for t in range(len(dis_t)):
        idx = type_ls.index(dis_t[t])
        min_dis = eval(min_ls[idx])
        max_dis = eval(max_ls[idx])
        if min_dis < dis[t] < max_dis:
            into_ele = [dis[t], dis_t[t], dis_idx[t]]
            useful_dis.append(into_ele)


    return useful_dis

#useful_dis = dis_func(bb,coor)

def bond_class(ls):  # Sum of all bonding energy
    ene_out = 0
    for q in ls:
        xb = q[0]
        typ = q[1]
        idx = typeb.index(typ)
        equilibrium = numpy.float64(xeqb[idx])
        x4 = numpy.float64(xconsb[0][idx])
        x3 = numpy.float64(xconsb[1][idx])
        x2 = numpy.float64(xconsb[2][idx])
        x1 = numpy.float64(xconsb[3][idx])
        cons = numpy.float64(xconsb[4][idx])
        dx = xb - equilibrium
        ener_ele = x4*dx**4 + x3*dx**3 + x2*dx**2 + x1*dx + cons
        ene_out = ene_out + ener_ele

    return ene_out


def hb_func(coor, useful_dis):
    cand_hdonor = []
    for j in range(0, len(useful_dis)):
        if useful_dis[j][1] == 'O-H':
            h_idx = eval(useful_dis[j][2].split('-')[1]) - 1
            cand_hdonor.append(h_idx)
        elif useful_dis[j][1] == 'H-O':
            hh_idx = eval(useful_dis[j][2].split('-')[0]) - 1
            cand_hdonor.append(hh_idx)
        else:
            continue

    cand_oacc = []
    for d in range(0, len(coor)):
        if coor[d][0] == 'O':
            cand_oacc.append(d)

    cand_oh = []
    for s in range(0, len(cand_hdonor)):
        h_coor = numpy.array(coor[cand_hdonor[s]][1:])
        for f in range(0, len(cand_oacc)):
            o_coor = numpy.array(coor[cand_oacc[f]][1:])
            dis_oh = numpy.linalg.norm(h_coor - o_coor)
            if 1.3 < dis_oh < 2.6:
                cand_oh.append(dis_oh)

    return cand_oh

#useful_hb = hb_func(coor, useful_dis)
def moorse_pot(ls):
    m_part = 0
    for s in ls: # Here I used Moorse potential to decribe the HB
        #lj_ele = (2.6/s)**12 - (2.6/s)**6
        m_ele = 0.659*(math.exp(-2*1.631*(s-2.347))-2*math.exp(-1.631*(s-2.347)))
        m_part = m_part + m_ele
    return m_part

def ang_func(bb, useful_dis, coor):
    #coor = Sugar(bb).next()
    useful_ang = []
    for d in range(len(useful_dis)):
        ii = useful_dis[d][2].split('-')
        for s in range(d+1, len(useful_dis)):
            jj = useful_dis[s][2].split('-')
            dup_idx = [x for x in ii if x in jj]
            if len(dup_idx) == 1:
                ii2 = ii.copy()
                ii2.remove(dup_idx[0])
                jj.remove(dup_idx[0])
                e1 = ii2[0]
                e2 = dup_idx[0]
                e3 = jj[0]
                e1e = eval(e1) -1
                e2e = eval(e2) -1
                e3e = eval(e3) -1
                e1c = numpy.array(coor[e1e][1:])
                e2c = numpy.array(coor[e2e][1:])
                e3c = numpy.array(coor[e3e][1:])
                ba = e1c - e2c
                bc = e3c - e2c
                cosine_angle = numpy.dot(ba, bc) / (numpy.linalg.norm(ba) * numpy.linalg.norm(bc))
                ang = numpy.arccos(cosine_angle)/3.14159*180
                if str(coor[e1e][0])=='Na' or str(coor[e2e][0])=='Na' or str(coor[e3e][0])=='Na':
                    continue
                else:
                    ang_a_type = str(coor[e1e][0]) + '-' + str(coor[e2e][0]) + '-' + str(coor[e3e][0])
                    idx_elea = str(e1) + '-' + str(e2) + '-' + str(e3)
                    into_elea = [ang, ang_a_type, idx_elea]
                    useful_ang.append(into_elea)
    return useful_ang

#useful_ang = ang_func(bb, useful_dis, coor)
#print(useful_ang)

def angle_class(ls): # Sum of all bonding energy
    ene_out = 0
    for q in ls:
        xa = q[0]
        typ = q[1]
        idx = typea.index(typ)
        equilibrium = numpy.float64(xeqa[idx])
        x4 = numpy.float64(xconsa[0][idx])
        x3 = numpy.float64(xconsa[1][idx])
        x2 = numpy.float64(xconsa[2][idx])
        x1 = numpy.float64(xconsa[3][idx])
        cons = numpy.float64(xconsa[4][idx])
        dx = xa - equilibrium
        ener_ele = x4*dx**4 + x3*dx**3 + x2*dx**2 + x1*dx + cons
        ene_out = ene_out + ener_ele
    return ene_out



def charges_func(bb, coor):
    ff = open(bb)
    gg = ff.readlines()
    filename = bb.split(".")[0]
    mull = single_cmd("grep -n 'Mulliken charges:' {}".format(bb))
    natom = single_cmd("grep -n 'NAtoms= ' {}".format(bb))
    num_atom = eval(natom[0].split()[2])

    starting_idx = eval(mull[0].split(":")[0]) + 1
    end_idx = eval(mull[0].split(":")[0]) + 1 + num_atom


    rows = []
    for i in range(starting_idx, end_idx):
        mull_ele = []
        mull_ele.append(gg[i])
        sub_ele = mull_ele[0].split()
        row_ele = sub_ele[2]
        rows.append(row_ele)
    final_ch_data = []
    for i in range(len(coor)):
        ai = numpy.array(coor[i][1:])
        ait = coor[i][0]
        chi = rows[i]
        for j in range(i+1, len(coor)):
            aj = numpy.array(coor[j][1:])
            cur_dis = numpy.linalg.norm(ai - aj)
            chj = rows[j]
            final_ch_data.append([cur_dis, chi, chj])

    return final_ch_data

#useful_charges = charges_func(bb, coor)
def coulomb_func(ls):
    coulomb_E = 0
    for d in range(len(ls)):
        c_E_ele = 8.988*pow(10, 19)*(eval(ls[d][1])*eval(ls[d][2]))/(ls[d][0])
        coulomb_E = coulomb_E + c_E_ele
    coulomb_E = coulomb_E*3.657*pow(10, -18)
    return coulomb_E

########################### Here is the general case ################
structures = file_list('.', keys='.out', recursive=False)
all_result = []
for u in structures:
    coor = Sugar(u).next()
    bb = basename(u)
    useful_dis = dis_func(bb, coor)
    useful_ang = ang_func(bb, useful_dis, coor)
    useful_hb = hb_func(coor, useful_dis)
    useful_charges = charges_func(bb, coor)
    BE = bond_class(useful_dis)
    AE = angle_class(useful_ang)
    HB = moorse_pot(useful_hb)
    CE = coulomb_func(useful_charges)
    sumup = BE + AE + HB + CE
    zpe = out2zpe(bb)
    #print([bb, BE, AE, HB, CE, sumup])
    all_result.append([bb, BE, AE, HB, CE, sumup, zpe])



print('--------------------')
print(all_result)
#print('Bonding E=',bond_class(useful_dis)) # This is all bonding energy in kcal/mol
#print('Angle E=',angle_class(useful_ang)) # This is all angle energy in kcal/mol
#print('Hydrogen Bonding=',moorse_pot(useful_hb)) # This is all Hydrogen bond described via Moorse potential
#print('Coulomb E=',coulomb_func(useful_charges)) # This is absolute value, only need to be recorded and compare, not added to the first 2
print('--------------------')


