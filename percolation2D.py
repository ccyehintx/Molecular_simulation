# Percolation simulation
# Author: Chieh-Chih (George) Yeh
# Date: 2022/10/29
# This script is created to study the phase transition of percolation theory
# in a 2D lattice
import numpy as np
import matplotlib.pyplot as plt
import random
import networkx as nx
N = 20 # Size of the matrix

size = N*N # Size of the multiplied matrix
zz = np.zeros((N, N)) # Creating a zero matrix for the multiplied matrix
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
#print(num_cluster(0.5, N)[0])
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
cls = num_cluster(0.5, N)[0]
bls = num_cluster(0.5, N)[3]
ori = show_bonding(cls, bls)[0]
des = show_bonding(cls, bls)[1]

xx = []
yy = []
probls = [0.1, 0.2, 0.3, 0.4, 0.5, 0.51, 0.53, 0.55, 0.57, 0.59, 0.6, 0.7, 0.8, 0.9]
for i in range(1, 10):
    prob = i*0.1
    nc = num_cluster(prob, N)
    nc1 = nc[1]
    nc2 = nc[2]
    x = range(nc1)#range(num_cluster(prob, N)[1])
    y = nc2#num_cluster(prob, N)[2]
    xx.append(prob)
    yy.append(nc1)
    #print(prob)
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    ax.set_ylim([0, size+1])
    ax.set_ylabel('Cluster')
    ax.set_xlabel('Number of nodes in each cluster')
    #fig, ax = plt.subplots()
    rects1 = ax.bar(x, y)
    title = 'Percolation with P={} and {} clusters'.format(prob, nc1)
    plt.title(title)
    fname = 'prob_{}.png'.format(i)
    plt.savefig(fname)
    plt.clf()
    print(num_cluster(prob, N)[1], num_cluster(prob, N)[2])

col_ls = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
#plt.xlabel('entry a')
#plt.ylabel('entry b')

ccount = 0
for i in range(len(ori)): # there's one more layer in this
    orix = []
    oriy = []
    desx = []
    desy = []
    if ccount > 5:
        colour = col_ls[-1]
    else:
        colour = col_ls[ccount]
    for j in range(len(ori[i])):
        orix.append(ori[i][j][0])
        oriy.append(ori[i][j][1])
        desx.append(des[i][j][0])
        desy.append(des[i][j][1])
    allx = orix + desx
    ally = oriy + desy
    #plt.scatter(allx, ally, color = colour)
    ccount = ccount + 1
#plt.show() ################# Uncomment to plot
plt.plot(xx, yy)
plt.title("Number of cluster v.s. Probability")
plt.xlabel("Probability")
plt.ylabel("Number of cluster")
plt.savefig('curve.png')


