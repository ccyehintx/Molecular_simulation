import numpy 
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from csv_tools import read_csv_file

k_r_dehy = read_csv_file("k_r_dehy.csv")["001_NUM"]
k_ro = read_csv_file("k_ro.csv")["001_NUM"]
k_l_dehy = read_csv_file("k_l_dehy.csv")["001_NUM"]
k_02A = read_csv_file("k_02A.csv")["001_NUM"]
k_03A = read_csv_file("k_03A.csv")["001_NUM"]

tspan = numpy.linspace(0, 0.01, 5000) ## time scale needs to be confirmed
yinit = [1, 0, 0, 0, 0]
k = [k_r_dehy[24], k_ro[24],  k_l_dehy[24], k_02A[24], k_03A[24]]

def f(t, y, k):
	dydt = [ -k[0]*y[0]-k[1]*y[0], k[0]*y[0]+k[2]*y[4], k[3]*y[4], k[4]*y[4], k[1]*y[0]-k[2]*y[4]-k[3]*y[4]-k[4]*y[4]]
	return dydt


sol = solve_ivp(lambda t, y: f(t, y, k),  [tspan[0], tspan[-1]], yinit, t_eval=tspan)


time = sol.t
y0 = sol.y[0]
y1 = sol.y[1]
y2 = sol.y[2]
y3 = sol.y[3]
y4 = sol.y[4]

plt.plot(time, y0, label='Ring glcA') 
plt.plot(time, y1, label='H2O') 
plt.plot(time, y2, label='02A') 
plt.plot(time, y3, label='03A') 
plt.plot(time, y4, label='Linear glcA')

#plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

plt.legend(bbox_to_anchor=(1, 1),
           bbox_transform=plt.gcf().transFigure)

plt.title('Concentration profile at 875K')
plt.xlabel('Time(s)')
plt.ylabel('Concentration ratio (-)')

plt.show()


