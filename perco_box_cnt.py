# Box counting method in a 2D percolation lattice
# Author: Chieh-Chih (George) Yeh
# Date: 2023/05/07
# This script is created to study the fractal character of percolation flow

import scipy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import networkx as nx
import math

def chk_region(idx,N):
    # This function tells the neighbours of the input index
    # cornors:
    ul = 0
    ur = N - 1
    ll = N*N - N
    lr = N*N - 1
    ulv = [1, N-1]
    urv = [N-1, 2*N-1]
    llv = [(N-1)*N - N, N*N - N + 1]
    lrv = [N*N - 2, (N-1)*N - 1]
    all_cor = [ul, ur, ll, lr]
    corv = [ulv, urv, llv, lrv]
    # upper 
    all_up = []
    all_down = []
    all_right = []
    all_left = []
    upv = []
    downv = []
    rightv = []
    leftv = []
    for i in range(1, N-1):
        all_up.append(i)
        upv.append([i-1, i+1, i+N])
        aright_idx = (i+1)*N-1
        all_right.append(aright_idx)
        rightv.append([aright_idx-1, aright_idx-N, aright_idx+N])
        adown_idx = N*N-N + i
        all_down.append(adown_idx)
        downv.append([adown_idx-1, adown_idx+1, adown_idx-N])
        aleft_idx = N*i
        all_left.append(aleft_idx)
        leftv.append([aleft_idx-N, aleft_idx+N, aleft_idx+1])
    all_spec = [all_cor, all_up, all_down, all_right, all_left]
    specv = [corv, upv, downv, rightv, leftv]
    for i in range(len(all_spec)):
        if idx in all_spec[i]:
            ridx = all_spec[i].index(idx)
            vv = specv[i][ridx]
            break
        else:
            vv = [idx - 1, idx + 1, idx - N, idx + N]
    return vv

def num_cluster(prob, N):
    bonding = []
    size = N*N
    bonding_id = []
    cluster_id = []
    cont = 0
    for i in range(size):
        nei = chk_region(i, N)
        group_id = -1
        lump = [i]
        for j in nei:
            rand = random.random()
            if rand < prob:
                if [j, i] not in bonding:
                    bonding.append([i, j])
                    #cont = cont + 1
                    lump.append(j)
                    bonding_id.append(cont)
                    cont = cont + 1

    for i in range(cont):
        group_id = -1
        if i == 0:
            cluster_id.append(bonding[i])
        else:
            for j in range(len(cluster_id)):
                if any(x in bonding[i] for x in cluster_id[j]):
                    group_id = j
                    break
            if group_id == -1:
                cluster_id.append(bonding[i])
            else:
                for k in bonding[i]:
                    cluster_id[group_id].append(k)

    v2cluster_id = []
    sep_count = []
    for i in cluster_id:
        ii = [*set(i)]
        sep_count.append(len(ii))
        v2cluster_id.append(ii)
    sep_count.sort()
    return v2cluster_id, len(v2cluster_id), sep_count, bonding

def show_bonding(cluster, bonding):
    ori = []
    des = []
    for i in cluster:
        ori.append(list())
        des.append(list())
        for k in i:
            j = bonding[k]
            j1r = j[0]//N
            j1c = (j[0]+1)%N - 1
            if j1c == -1:
                j1c = j1c + N
            j1 = [j1r, j1c]
            j2r = j[1]//N
            j2c = (j[1]+1)%N - 1
            if j2c == -1:
                j2c = j2c + N
            j2 = [j2r, j2c]
            ori[-1].append(j1)
            des[-1].append(j2)
    return ori, des


def bd_nodes(N, w):
    # This function will give the boundary nodes for breakdown
    dim = w + 1
    zz = np.zeros((dim, dim))
    for i in range(w+1):
        ii = i*((N-1)/w)*N
        for j in range(w+1):
            jj = j*((N-1)/w)
            ij = ii + jj
            zz[i][j] = ij
    return zz
    
def int_gen(ls, N):
    # This provide the list of indices include in the breakdown box
    ul = ls[0]
    ur = ls[1]
    ll = ls[2]
    lr = ls[3]
    r1 = int(ul)%N
    r2 = math.ceil(ur)%N
    row_ls = range(r1, r2+1)
    col_ls = []
    col_ls.append(int(ul/N)*N)
    uul = col_ls[0] + N
    lll = int(ll/N)*N
    while uul < lll:
        col_ls.append(uul)
        uul = uul + N
    col_ls.append(lll)
    int_ls = []
    for i in col_ls:
        for j in row_ls:
            ij = i + j
            int_ls.append(ij)
    return int_ls

# create box index from the box_node matrix
def box_cnt(bd, lss, N, w):
    # This function counts how many boxes needed 
    nodes = lss.copy()
    box_count = 0
    see = []
    for i in range(w):
        for j in range(w):
            ul = bd[i][j]
            ur = bd[i][j+1]
            ll = bd[i+1][j]
            lr = bd[i+1][j+1]
            ls = [ul, ur, ll, lr]
            inttt = int_gen(ls, N)
            common_list = list(set(inttt).intersection(nodes))
            if len(common_list) > 0:
                box_count = box_count + 1
                for k in common_list:
                    nodes.remove(k)
    return box_count

# Models to be fitted
def power_law(l, a, d):
    return a*l**(-d)

def poisson_dis(l, a, d):
    #ff = a*math.exp(-l)*l**d/((2*3.14159*d)**0.5*(d/math.exp(1))**d)
    ff = a*math.exp(-d)*d**l/((2*3.14159*l)**0.5*(l/math.exp(1))**l)
    return ff

# Getting R square
def rsquare(number_of_boxes, popt, l, dis):
    for k in range(len(l)):
        residuals = []
        rr = number_of_boxes[k] - dis(l[k],popt[0],popt[1])
        residuals.append(rr**2)
    squaresumofresiduals = np.sum(residuals)
    squaresum = np.sum((number_of_boxes-np.mean(number_of_boxes))**2)
    R2 = 1 - (squaresumofresiduals/squaresum)
    return R2

############ Input Parameters ############
N = 40 # Size of the matrix
wls = [2,30] # Split the boundary into w pieces
prob = 0.5 # Probability of the flow

cls = num_cluster(prob, N)[0]

concat_list = [j for i in cls for j in i]
lss = list(set(concat_list))
number_of_boxes = []
l = []
for w in range(wls[0], wls[1]):
    bd = bd_nodes(N,w)
    bc = box_cnt(bd, lss, N, w)
    number_of_boxes.append(bc)
    l.append(N/w)
    #print('The box is divided into', w)
    #print('The box length is then', N/w)
    #print('The number of boxes required', bc)

poptpl, pcovpl = scipy.optimize.curve_fit(f=power_law, xdata=l, ydata=number_of_boxes, p0=(30, 2))

poptpd, pcovpd = scipy.optimize.curve_fit(f=poisson_dis, xdata=l, ydata=number_of_boxes, p0=(30, 2))

pl_y = []
for i in l:
    pl_y.append(power_law(i, poptpl[0], poptpl[1]))
pd_y = []
for i in l:
    pd_y.append(poisson_dis(i, poptpd[0], poptpd[1]))

R2pl = rsquare(number_of_boxes, poptpl, l, power_law)
R2pd = rsquare(number_of_boxes, poptpl, l, poisson_dis)

plt.scatter(l, number_of_boxes, label='Percolation simulation')
plt.plot(l, pl_y, color='red', label='Power law')
plt.plot(l, pd_y, color='green', label='Poisson distribution')
plt.xlabel('Box Length')
plt.ylabel('Number of Boxes')
plt.legend()
plt.show()


print('----------Computed results-------------')
print('The 2D lattice with the boundary length of {}'.format(N))
print('The probability of the percolation is defined as {}'.format(prob))
print('The lattice is divided into pieces from {} to {}'.format(wls[0],wls[1]))
print('Power law fitted fractal dimension',poptpl[1])
print('Poisson distribution fitted fractal dimension',poptpd[1])
print('Power law R squared value', R2pl)
print('Poisson distribution R squared vaule', R2pd)
