import numpy as np
import math
import sys
import scipy.constants as scc
wham_file = sys.argv[1]

#kT = 0.6 
kT = 0.5921868 
bin_centers = []
free_energy = []
fe_err = []
C_values = []
with open(wham_file,'r') as f:
        for line in f:
                temp = line.split()
                if temp[0] == '#Coor' or temp[0] == '#Window':
                        continue
                elif temp[1] == 'inf' or temp[2] == 'nan':
                        continue
                elif temp[0][0] != '#':
                        bin_centers.append(float(temp[0]))
                        free_energy.append(float(temp[1]))
                        fe_err.append(float(temp[2]))
                # GRABBING DATA TO BE USED TO CREATE THE UNSTICKED FE PLOTS; THIS ELSE STATEMENT GRABS THE CONSTANT VALUES THAT IS USED TO ALIGN THE WINDOWS WITH EACH OTHER;                                       
                else:
                        C_values.append(float(temp[1]))

bin_centers = np.array(bin_centers)
free_energy = np.array(free_energy)
fe_err = np.array(fe_err)
C_values = np.array(C_values)

min_index = np.argmin(free_energy)
max_index = len(bin_centers) - 1
PMF_max = free_energy[max_index]
print "Gw = %s kCal/mol, %s kJ/mol" %(-PMF_max-1.2*math.log(bin_centers[max_index]/bin_centers[min_index]),(-PMF_max-1.2*math.log(bin_centers[max_index]/bin_centers[min_index]))*4.184)
print fe_err[2]
R_max = bin_centers[max_index]

delta_PMF = -kT*np.log(4*np.pi*R_max*R_max) - PMF_max
#print "delta_PMF = %s" %(delta_PMF)

#print free_energy
free_energy += delta_PMF

#print free_energy
exp_free_energy = np.exp(-free_energy/kT)
#for i in range(len(bin_centers)):
#	print bin_centers[i], exp_free_energy[i]
Area = 0
delta_x = (bin_centers[len(bin_centers)-1] - bin_centers[0])/len(bin_centers)
Area_error=0
#print bin_centers[max_index], free_energy[max_index]
#print delta_x
for i in range(len(bin_centers)):
    if i == 0 or i == len(bin_centers)-1:
        Area +=  exp_free_energy[i]
	Area_error += exp_free_energy[i]*exp_free_energy[i]*fe_err[i]*fe_err[i]/kT/kT

    elif i%2 == 0:
        Area += 2*exp_free_energy[i]
	Area_error += 2*2*exp_free_energy[i]*exp_free_energy[i]*fe_err[i]*fe_err[i]/kT/kT
    else:
        Area += 4*exp_free_energy[i]
	Area_error += 4*4*exp_free_energy[i]*exp_free_energy[i]*fe_err[i]*fe_err[i]/kT/kT

Area *= (delta_x)/3
Area_error = math.sqrt(Area_error)
Area_error *= (delta_x)/3
C = scc.N_A*10**-27
Kd = Area*C
Kd_error = Area_error*C

print "Kd = %e" %(Kd)
print "Kd error = %e" %(Kd_error)
Gd = -kT*np.log(Kd)*4.184
Gd_error = kT*4.184*Kd_error/Kd
print "Gd = %s" %(Gd)
print "Gd error = %s" %(Gd_error)
