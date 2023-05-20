# This script will turn the time series into singular spectrum

import matplotlib.pyplot as plt
import numpy as np
        
def f_2sq(s, prof):
    l = len(prof)
    ls = int(l/s)
    v_ls = []
    for v in range(ls):
        i0 = int(v*s)
        iff = int((v+1)*s)
        cur_ls = prof[i0:iff]
        x = list(range(len(cur_ls)))
        y_curl = np.polyfit(x, cur_ls, 3) # give 4 values
        y1, y2, y3, y4 = y_curl
        cur_v = 0
        for i in range(len(cur_ls)):
            ii = (cur_ls[i] - y1*(i**3) - y2*(i**2) - y3*i - y4)**2
            cur_v = cur_v + ii
        v_ls.append(cur_v/s)
    return v_ls
     
def f_q(f2, s, q):
    ls = len(f2)
    allf = 0
    for f in f2:
        ff = f**(q/2)
        allf = allf + ff
    result = (1/ls*allf)**(1/q)
    return result

def hs_con(fq_ls, s_ls):
    fq_ln = []
    s_ln = []
    for i in range(len(s_ls)):
        fq_ln.append(np.log(fq_ls[i]))
        s_ln.append(np.log(s_ls[i]))
    gr, con = np.polyfit(s_ln, fq_ln, 1)
    return gr

def tau_ls(hs, q):
    return hs*q - 1

def alpha_fd(x, y):
    # apply finite difference method
    # y is tau list and x is q list
    alpha_ls = []
    for i in range(len(x)-1):
        x0, x1 = x[i], x[i+1]
        y0, y1 = y[i], y[i+1]
        alpha = (y1 - y0)/(x1 - x0)
        alpha_ls.append(alpha)
        # grab all q except the first one
    return alpha_ls

def extract_file(file, si, fi, ex_id):
    # 152, 1152, 2
    sf = fi - si
    mylines = []
    with open(file, 'rt') as myfile:
        for myline in myfile:
            mylines.append(myline)
    x = range(0, sf)
    real_num = []
    for i in range(si, fi):
        iv = mylines[i].split()[ex_id]
        ii = round(eval(iv), 4)
        real_num.append(ii)
    return real_num

def profile_ls(real_num):
    avg = sum(real_num)/len(real_num)
    profile = []
    cur_sum = 0
    for i in real_num:
        iele = i - avg + cur_sum
        profile.append(iele)
        cur_sum = iele
    return profile

raw_data = extract_file('log.lammps', 152, 1152, 2)
profile = profile_ls(raw_data)

q_ls = list(np.arange(-2, 4, 0.3))
s_ls = [5, 10, 20, 40, 50, 100, 200, 250]
all_tau = []
for q in q_ls:
    cnt_fq = []
    for s in s_ls:
        f2 = f_2sq(s, profile)
        fq = f_q(f2, s, q)
        cnt_fq.append(fq)
    hs = hs_con(cnt_fq, s_ls)
    tau = tau_ls(hs, q)
    all_tau.append(tau)
al_ls = alpha_fd(q_ls, all_tau)

f_ls = []
for a in range(len(al_ls)):
    aa = a + 1
    ff = al_ls[a]*q_ls[aa] - all_tau[aa]
    f_ls.append(ff)
    #print(ff, al_ls[a], q_ls[aa])

plt.plot(al_ls, f_ls, label='Original')
plt.xlabel('Singularity Strength')
plt.ylabel('Singularity Spectrum')
plt.title('Multifractal Spectrum')
#plt.legend()
plt.show()
#plt.savefig('q1.png')
