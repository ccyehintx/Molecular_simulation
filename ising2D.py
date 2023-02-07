# This is a script to perform 2D ising model
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import moviepy.video.io.ImageSequenceClip
N = 50
timeN = 500
sites = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        k = 2*random.randint(0, 1) - 1
        sites[i, j] = k
def metMC(sites, N, b):
    sitesz = sites.copy()
    for i in range(N):
        for j in range(N):
            x = random.randint(0, N-1)
            y = random.randint(0, N-1)
            s = sitesz[x, y]
            xp, xm = (x+1)%N, (x-1)%N # to accout for periodicity
            yp, ym = (y+1)%N, (y-1)%N
            ss = sitesz[xp, y] + sitesz[xm, y] + sitesz[x, yp] + sitesz[x, ym]
            e = 2*s*ss
            if e < 0:
                s = s*(-1)
            elif random.random() < np.exp(-e*b):
                s = s*(-1)
            sitesz[x, y] = s
    return sitesz
#b = 1/3.28
b = 1/1.5
M = []
M.append(sites)
for i in range(timeN):
    mm = metMC(M[-1], N, b)
    M.append(mm)
fname = []
for i in range(timeN):
    if i%5 == 0:
        mat = M[i]
        plt.matshow(mat, cmap=plt.cm.Blues)
        plt.title("time = {}".format(i))
        plt.savefig('num_{}.png'.format(i))
        fname.append('num_{}.png'.format(i))
#plt.matshow(M[-1], cmap=plt.cm.Blues)
#plt.show()
image_folder='./'
fps=1

image_files = [os.path.join(image_folder,img)
        for img in fname]
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
clip.write_videofile('my_video.mp4')
