# Faster way to calculate coordinates by matrix operations

import numpy as np
import gzip

linker_num = 3000
colloid_num = 1000

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

nstep = 100
sidx = 31009*nstep
fidx = sidx + 31008
sort_info = []
for i in range(sidx+9, fidx+1):
    sort_info.append(raw_line(mylines[i]))

patch_ls = []
head_ls = []
for i in range(len(sort_info)):
    if sort_info[i][1] == 2:
        patch_ls.append(i) # this append the index of sort_info
    elif sort_info[i][1] == 3:
        head_ls.append(i)


bond_ls = []
for b in range(linker_num):
    bond_ls.append(list())
head_all_info = []
head_iidx = []
for j in head_ls:
    head_all_info.append(sort_info[j][2])
    head_iidx.append(sort_info[j][0])
#head_coor = np.array(head_all_info[:,2])
for i in patch_ls:
    cidx = (sort_info[i][0] - 1)//7
    icoor = sort_info[i][2]
    collect_dis = dis_mat(icoor, head_all_info)
    test_list = list(collect_dis)
    res = [idx for idx, val in enumerate(test_list) if val < 0.5]
    if len(res) > 0:
        for k in res:
            bbidx = head_iidx[k]
            jjidx = (bbidx - 7000 - 1)//8
            bond_ls[jjidx].append(cidx)

#print(bond_ls)
cnt = 0
for i in bond_ls:
    cnt = cnt + len(i)
        

print(cnt)
