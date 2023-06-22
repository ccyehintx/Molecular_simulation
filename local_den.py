# IPython log file

import subprocess
import numpy as np 

def raw_line(line):
    info = []
    info.append(eval(line.split(' ')[0]))
    info.append(eval(line.split(' ')[1]))
    coor = np.array([eval(i) for i in line.split(' ')[2:5]])
    info.append(coor)
    return info


def ide_box(brange, unitlen, coor):
    nsubs = (abs(brange[1] - brange[0])/unitlen)
    tns = nsubs**3
    x,y,z = coor - np.array([0.001, 0.001, 0.001])
    zbox = (z - brange[0])//unitlen
    ybox = (y - brange[0])//unitlen
    xbox = (x - brange[0])//unitlen
    idx = xbox + ybox*nsubs + zbox*nsubs**2
    return idx


brange = [0, 20]
unitlen = 2
tns = int((abs(brange[1] - brange[0])/unitlen)**3)
#ttls = np.zeros(tns)
filename = 'argon.lammpstrj'
mylines = []
with open(filename, mode ='r')as file:
    for ll in file:
        mylines.append(ll)

cmd = "grep -n 'ITEM: NUMBER OF ATOMS' {}".format(filename)
command = subprocess.check_output(cmd, shell=True)
cc = command.decode("utf-8").split("\n")
for i in cc:
    if len(i) > 3:
        nline = eval(i.split(':')[0])
        natom = eval(mylines[nline])
        sidx = nline + 6
        eidx = sidx + natom
        ttlist = np.zeros(tns) # This list gives the num of colloids in each small box
        for j in range(sidx, eidx):
            coor = raw_line(mylines[j])[2]
            boxidx = int(ide_box(brange, unitlen, coor))
            if boxidx > len(ttlist):
                print(coor, boxidx)
            else:
                ttlist[boxidx] = ttlist[boxidx] + 1
        std = np.std(ttlist)
        print(std)

