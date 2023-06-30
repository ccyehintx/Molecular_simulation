# IPython log file

import subprocess
import numpy as np
import gzip
import matplotlib.pyplot as plt
import random
import networkx as nx
import matplotlib.cm as cm
from matplotlib.colors import Normalize
mylines = []
with gzip.open('quench.data.gz', mode ='r')as file:
    for ll in file:
        mylines.append(ll)
        
def dis_mat(coor_patch, coor_heads):
    A = np.ones((len(coor_heads), 3))
    A[:, 0] = A[:, 0]*coor_patch[0]
    A[:, 1] = A[:, 1]*coor_patch[1]
    A[:, 2] = A[:, 2]*coor_patch[2]
    B = coor_heads
    C = B - A
    D = np.square(C)
    E = np.sum(D, axis=1)
    F = np.sqrt(E)
    return F
def raw_line(line):
    info = []
    info.append(eval(line.split()[0]))
    info.append(eval(line.split()[2]))
    coor = np.array([eval(i) for i in line.split()[3:6]])
    info.append(coor)
    return info

def sort_info_func(mylines, sidx, fidx):
    sort_info = []
    for i in range(sidx+9, fidx+1):
        sort_info.append(raw_line(mylines[i]))
    return sort_info

def pandh_ls(sort_info):
    patch_ls = []
    head_ls = []
    for i in range(len(sort_info)):
        if sort_info[i][1] == 2: # 2 is patch
            patch_ls.append(i)
        elif sort_info[i][1] == 3: # 3 is head
            head_ls.append(i)
    return patch_ls, head_ls
def bond_ls_func(sort_info, linker_num, patch_ls, head_ls):
    bond_ls = []
    for b in range(linker_num):
        bond_ls.append(list())
    head_all_info = []
    head_iidx = []
    for j in head_ls:
        head_all_info.append(sort_info[j][2])
        head_iidx.append(sort_info[j][0])
    for i in patch_ls:
        cidx = (sort_info[i][0] - 1)//7
        icoor = sort_info[i][2]
        collect_dis = dis_mat(icoor, head_all_info)
        test_list = list(collect_dis)
        res = [idx for idx, val in enumerate(test_list) if val < 0.5]
        if len(res) > 0:
            for k in res:
                bbidx = head_iidx[k]
                jjidx = (bbidx - 7000 - 1)//8 #(1000 + 6*1000 = 7000)
                bond_ls[jjidx].append(cidx)
    return bond_ls
def gen_G(size, bonding):
    # Generate a graph G
    G = nx.Graph()
    for i in range(size):
        G.add_node(i)
    for i in bonding:
        G.add_edge(i[0],i[1])
    return G

def avg_deg(degree_ls):
    # Generate average degree
    sum_deg = []
    for i in degree_ls:
        sum_deg.append(i[1])
    return sum(sum_deg)/len(sum_deg)

def global_eff(G, bigc, size):
    p = nx.shortest_path(G)
    node_idx = -1
    for i in range(size):
        not_in_the_same_c = []
        for j in range(len(bigc)):
            if i in bigc[j]:
                node_idx = j
    return p
def avg_nodal_conn(G, size):
    sum_k = []
    for i in range(size):
        for j in range(i+1,size):
            sum_k.append(nx.node_connectivity(G, i, j))
    avg_nc = 2*sum(sum_k)/size/(size-1)
    return avg_nc

def locate_coor(idx, n):
    z = idx//(n**2)
    ii = idx%(n**2)
    x = ii%n
    y = ii//n
    return np.array([x, y, z])

def ide_box(brange, unitlen, coor):
    nsubs = (abs(brange[1] - brange[0])/unitlen)
    tns = nsubs**3
    x,y,z = coor - np.array([0.001, 0.001, 0.001])
    zbox = (z - brange[0])//unitlen
    ybox = (y - brange[0])//unitlen
    xbox = (x - brange[0])//unitlen
    idx = xbox + ybox*nsubs + zbox*nsubs**2
    return idx

brange = [-55, 55]
unitlen = 11
tns = int((abs(brange[1] - brange[0])/unitlen)**3)
n = int((abs(brange[1] - brange[0])/unitlen))
linker_num = 3000
colloid_num = 1000
npatch = 6
nbead = 6

#n = 10 ############
graph_collect = []
for tstep0 in range(10, 220, 10):
    #tstep = 160
    tstep = int(tstep0)
    sidx = 31009*tstep
    fidx = sidx + 31008
    sort_info = sort_info_func(mylines, sidx, fidx)
    patch_ls, head_ls = pandh_ls(sort_info)
    bond_ls = bond_ls_func(sort_info, linker_num, patch_ls, head_ls)
    realbond = []
    for i in bond_ls:
        if len(i) == 2:
            realbond.append(i)        
    G = gen_G(colloid_num, realbond)
    close_cen = nx.closeness_centrality(G)
    between_cen = nx.betweenness_centrality(G, k=None, normalized=True, weight=None, endpoints=False, seed=None)
    
    bblist = list(between_cen.values())
    cclist = list(close_cen.values())
    tblist = np.zeros(tns) # this is the list of betweenness cen
    tclist = np.zeros(tns) # this is the list of closeness cen
    tdlist = np.zeros(tns)
    for j in sort_info:
        if j[1] == 1: #this is colloid
            cidx = j[0]//7
            coor = j[2]
            boxidx = int(ide_box(brange, unitlen, coor))
            if boxidx > len(tblist): # this is calling an error
                print(coor, boxidx)
            else:
                tblist[boxidx] = tblist[boxidx] + bblist[cidx]
                tclist[boxidx] = tclist[boxidx] + cclist[cidx]
                tdlist[boxidx] = tdlist[boxidx] + 1
           
    cmb = np.zeros((n, n))
    cmc = np.zeros((n, n))
    cmd = np.zeros((n, n))
    for i in range(len(tblist)):
        x, y, z = locate_coor(i, n)
        cmb[y][x] = cmb[y][x] + tblist[i]
        cmc[y][x] = cmc[y][x] + tclist[i]
        cmd[y][x] = cmd[y][x] + tdlist[i]
    edge_list = G.edges()
    degree_ls = nx.degree(G)
    big_cluster = list(nx.connected_components(G))
    bigc_len = []
    for i in big_cluster:
        bigc_len.append(len(i))
    ge = nx.global_efficiency(G)
    cc = nx.clustering(G)
    sum_cc = []
    for i in range(len(cc)):
        sum_cc.append(cc[i])
    #assort_coe = nx.degree_assortativity_coefficient(G)
    allres = [tstep, 2*len(edge_list)/colloid_num/(colloid_num-1), avg_deg(degree_ls), max(bigc_len), ge, sum(sum_cc)/len(sum_cc)]
    graph_collect.append(allres)
    print(allres)
    #plt.imshow(cmb, vmin=0, vmax=0.5, interpolation='none')
    #plt.colorbar()
    #plt.title('Betweenness Centrality at {}'.format(tstep))
    #plt.savefig('bc_{}.png'.format(tstep))
    #plt.clf()
    #plt.imshow(cmd, interpolation='none')
    #plt.colorbar()
    #plt.title('Local Density at {}'.format(tstep))
    #plt.title('Closeness Centrality at {}'.format(tstep))
    #plt.savefig('ld_{}.png'.format(tstep))
    #plt.clf()
#print(graph_collect)
