from os_tools import single_cmd
from os_tools import file_list
import numpy
import pandas as pd
import math
from gaussian_tools import out2zpe
from csv_tools import write_csv_file

ts_file = file_list(".", keys=".out TS", recursive=False)
global_min = file_list(".", keys=".out -TS", recursive=False)
zpe = []
k_master = []
Temperature = []

##all necessary inputs

kb = 1.38064852e-23
h = 6.62607004e-34
Ti = 298.15
Temperature.append(Ti)
for e in range(300, 2025, 25):
	Temperature.append(e)
glo_zpe = out2zpe(global_min[0])
sym_key = single_cmd("grep -n 'Rotational symmetry' {}".format(global_min[0]))
sym_no = int(sym_key[0].split(".\n")[0].split(" ")[-1])

mass_key = single_cmd("grep -n 'Molecular' {}".format(global_min[0]))
mole_mass = eval(mass_key[0].split("amu")[0].split(" ")[-2])
amu_to_kg = 1.66053904e-27
m = mole_mass*amu_to_kg

#########Here will add another temperature loop
for t in Temperature:
	#calculating the partition function of global_min first
	kx_this_temp = [] ##will sum up all k in this temp later
	fre_key_g = single_cmd("grep -n 'Frequencies' {}".format(global_min[0]))
	f_g = open(global_min[0])
	dd_g = f_g.readlines()
	final_fre_g = []
	for e in fre_key_g:
    		num = e.split("--")[1].split("\n")[0].split()
    		for y in num:
        		final_fre_g.append(y)

	vib_ele_g = []
	for j in final_fre_g:
    		s = 1/(1-math.exp(-h*eval(j)*2.998e10/kb/t))
    		vib_ele_g.append(s)

	qvib_g_first = numpy.prod(vib_ele_g[0:20])
	qvib_g_second = numpy.prod(vib_ele_g[20:40])
	qvib_g_third = numpy.prod(vib_ele_g[40:])

	moment_key_g = single_cmd("grep -n 'Principal axes' {}".format(global_min[0]))
	moment_info_g = dd_g[int(moment_key_g[0].split(":")[0])+1]
	gg_g = moment_info_g.split("\n")[0].split("--  ")[1]
	ix_g = eval(gg_g.split(".")[0]+"."+gg_g.split(".")[1][0:5])
	iy_g = eval(gg_g.split(".")[1][5:]+"."+gg_g.split(".")[2][0:5])
	iz_g = eval(gg_g.split(".")[2][5:]+"."+gg_g.split(".")[3])
	rx_g = h**2/(8*3.14159**2*ix_g*kb)
	ry_g = h**2/(8*3.14159**2*iy_g*kb)
	rz_g = h**2/(8*3.14159**2*iz_g*kb)

	qr_g = 3.14159**0.5/sym_no*(t**(3/2)/((rx_g*ry_g*rz_g)**0.5))


	#calculating for the partition functions of ts
	for p in ts_file:
		f = open(p)
		dd = f.readlines()
		ts_zpe = out2zpe(p)
		fre_key = single_cmd("grep -n 'Frequencies' {}".format(p))

		final_fre = []
		for w in fre_key:
    			num = w.split("--")[1].split("\n")[0].split()
    			for y in num:
        			final_fre.append(y)

		del final_fre[0]
		vib_ele = []
		for f in final_fre:
    			s = 1/(1-math.exp(-h*eval(f)*2.998e10/kb/t))
    			vib_ele.append(s)
		qvib_first = numpy.prod(vib_ele[0:20])
		qvib_second = numpy.prod(vib_ele[20:40])
		qvib_third = numpy.prod(vib_ele[40:])

		moment_key = single_cmd("grep -n 'Principal axes' {}".format(p))
		moment_info = dd[int(moment_key[0].split(":")[0])+1]
		gg = moment_info.split("\n")[0].split("--  ")[1]
		ix = eval(gg.split(".")[0]+"."+gg.split(".")[1][0:5])
		iy = eval(gg.split(".")[1][5:]+"."+gg.split(".")[2][0:5])
		iz = eval(gg.split(".")[2][5:]+"."+gg.split(".")[3])
		rx = h**2/(8*3.14159**2*ix*kb)
		ry = h**2/(8*3.14159**2*iy*kb)
		rz = h**2/(8*3.14159**2*iz*kb)

		qr = 3.14159**0.5/sym_no*(t**(3/2)/((rx*ry*rz)**0.5))
		qv_ratio_1 = qvib_first/qvib_g_first
		qv_ratio_2 = qvib_second/qvib_g_second
		qv_ratio_3 = qvib_third/qvib_g_third
		qr_ratio = qr/qr_g
		q_all_ratio = qv_ratio_1*qv_ratio_2*qv_ratio_3*qr_ratio

		Eax = (ts_zpe - glo_zpe)*2625.5
		kx = (kb*t/h)*(q_all_ratio)*math.exp(-Eax*1.6605e-21/kb/t)
		kx_this_temp.append(kx)

	kxx = sum(kx_this_temp)
	k_master.append(kxx)

rows = []
for e in range(0, len(k_master)):
	rows.append([Temperature[e], " ", k_master[e]])
cols = ["Temperature(K)", " ", "Cumulative rate constant(1/s)"]
data = pd.DataFrame(columns=cols)
data[cols] = pd.DataFrame(rows)
output_name = "Rate_constant_profile.csv"
write_csv_file(output_name, data)

