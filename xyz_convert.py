# IPython log file
import matplotlib.pyplot as plt
import pandas as pd
import math
import cmath
import numpy as np
import networkx as nx
def raw_line(line):
    coor = np.array([eval(i) for i in line.split(' ')[1:4]])
    return coor
def sort_info_func(mylines, sidx, fidx):
    an = int(fidx - sidx)
    sort_info = np.zeros((an, 3))
    for j in range(sidx, fidx):
        jj = j - sidx
        sort_info[jj][0] = raw_line(mylines[j])[0]
        sort_info[jj][1] = raw_line(mylines[j])[1]
        sort_info[jj][2] = raw_line(mylines[j])[2]
    return sort_info

def dis_mat(coor0, coorls):
    A = np.ones((len(coorls), 3))
    A[:, 0] = A[:, 0]*coor0[0]
    A[:, 1] = A[:, 1]*coor0[1]
    A[:, 2] = A[:, 2]*coor0[2]
    B = coorls
    C = B - A
    D = np.square(C)
    E = np.sum(D, axis=1)
    F = np.sqrt(E)
    return F

def gen_G(size, bonding):
    # Generate a graph G
    G = nx.Graph()
    for i in range(size):
        G.add_node(i)
    for i in bonding:
        G.add_edge(i[0],i[1])
    return G
def collect_info(atom_num, tr):
    sidx = 2 + tr*(2+atom_num)
    fidx = sidx + atom_num
    sort_info = sort_info_func(mylines, sidx, fidx)
    bond_ls = []
    for i in range(len(sort_info)):
        ff = dis_mat(sort_info[i], sort_info)
        test_list = list(ff)
        res = [idx for idx, val in enumerate(test_list) if val < 2.2 and val > 0.001]
        if len(res) > 0:
            for cc in res:
                if [cc, i] not in bond_ls:
                    bond_ls.append([i, cc])
    G = gen_G(atom_num, bond_ls)
    return G, sort_info

def psi_calc(G, sort_info):
    collect_psi = []
    atom_num = len(sort_info)
    for i in range(atom_num):
        nls = [ele for ele in G.neighbors(i)]
        if len(nls) > 0:
            A = np.ones((len(nls), 3))
            A[:, 0] = A[:, 0]*sort_info[i][0]
            A[:, 1] = A[:, 1]*sort_info[i][1]
            A[:, 2] = A[:, 2]*sort_info[i][2]
            B = np.zeros((len(nls), 3))
            for ii in nls:
                iidx = nls.index(ii)
                B[iidx][0] = sort_info[iidx][0]
                B[iidx][1] = sort_info[iidx][1]
                B[iidx][2] = sort_info[iidx][2]
            r_ij_m = B - A
            psi = 0
            for rm in r_ij_m:
                tij = np.arctan2([rm[1]], [rm[0]])
                cx = 6*tij
                ccx = complex(0, cx)
                psi = psi + cmath.exp(ccx)/6
            ppsi = ((psi.real)**2 + (psi.imag)**2)**0.5
            collect_psi.append(ppsi)
            #collect_psi.append(psi/len(nls))
        else:
            collect_psi.append(0)
    return collect_psi

filename = 'mag10_eq.xyz'
mylines = []
with open(filename, mode ='r')as file:
    for ll in file:
        mylines.append(ll)
atom_num = eval(mylines[0])
traj_num = int(len(mylines)/(eval(mylines[0]) + 2))

tr = 200
G = collect_info(atom_num, tr)[0]
sort_info = collect_info(atom_num, tr)[1]

big_cluster = list(nx.connected_components(G))
bigc_len = []
for i in big_cluster:
    bigc_len.append(len(i))
def pdf_calc(bigc_len, nbin):
    count, bins_count = np.histogram(bigc_len, bins=nbin)
    pdf = count#/sum(count)
    xb = []
    yb = []
    for i in range(len(pdf)):
        if pdf[i] > 0:
            yb.append(pdf[i])
            xb.append(bins_count[1:][i])
    return xb, yb
xb, yb = pdf_calc(bigc_len, 20)
#plt.scatter(xb, yb, color="red", label="PDF")
plt.plot(xb, yb, '-o', color="red", label="PDF")
plt.legend()
plt.show()
#collect_psi = psi_calc(G, sort_info)
#print(collect_psi)
#print(sum(collect_psi)/len(collect_psi))
