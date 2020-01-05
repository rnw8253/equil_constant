## USAGE:
## python ka.py [pmf_file] [r1_value] [r2_value]
## NOTE: set r2_value = -1 to select maximum R distance for given concentration.

import numpy as np
import math
import sys
import scipy.constants as scc
from bisect import bisect_left


class Equil():

    def readWHAM(self, wham_file):
        bin_centers = []
        free_energy = []
        fe_err = []
        C_values = []
        with open(wham_file,'r') as f:
            for line in f:
                temp = line.split()## split each element of line into its own element of a list
                if temp[0] == '#Coor' or temp[0] == '#Window':
                    continue
                    ## just skip to next line if the line starts with #Coor or #Window, ie. the start of the pmf and pmf_vert listings
                elif temp[1] == 'inf' or temp[2] == 'nan':
                    continue
                    ## if the free energy of a line is infinite (ie. the second column = 'inf'), skip that line as well
                elif temp[0][0] != '#':## if the first character of the first element of a line is not '#' then do the following
                    bin_centers.append(float(temp[0]))## add the first column (coordinate value) to this array
                    free_energy.append(float(temp[1]))## add the second column (Free energy) to this array
                    fe_err.append(float(temp[2]))## add the third column (+/- in FE) to this array
                # GRABBING DATA TO BE USED TO CREATE THE UNSTICKED FE PLOTS; THIS ELSE STATEMENT GRABS THE CONSTANT VALUES THAT ARE USED TO ALIGN THE WINDOWS WITH EACH OTHER;                                       
                else:## the only values left are the vertical offset values at the end of the pmf file
                    C_values.append(float(temp[1]))## grap the second column (FE offset) and add it to this array

        ###### Turn these lists into arrays
        bin_centers = np.array(bin_centers)
        free_energy = np.array(free_energy)
        delta_x = bin_centers[1] - bin_centers[0] ## step size in Angstroms

        fe_err = np.array(fe_err)
        C_values = np.array(C_values)

        return bin_centers, free_energy, delta_x;



    def volumeCorrect_addDistance(self, bin_centers, free_energy, delta_x, max_dist):

        delta_PMF = -kT*np.log(bin_centers[-1]*bin_centers[-1]) - free_energy[-1] ## Why the 4*pi? 4*pi should be in every baseline correction but it only ends up shifting the entire PMF curve by log(4*pi) and since the FE scale is relative it doesn't matter here. Baseline entire pmf so that PMF[-1] is at: FE = 0.
        
        ## VOLUME CORRECT FREE ENERGY
        for i in range(len(free_energy)):
            free_energy[i] += 2.0*kT*np.log(bin_centers[i])

        free_energy += delta_PMF ## vertically shift every FE value of the PMF by the volume corrected value of the PMF at long R, (delta_PMF).

        # add distance points till it gets to max_dist
        while True:
            bin_centers = np.append(bin_centers, bin_centers[-1]+delta_x)
            free_energy = np.append(free_energy, 0.)
            if bin_centers[-1] >= max_dist:
                break

        exp_free_energy = np.exp(-free_energy/kT)## make Boltzmann factor for each free_energy. g(r) = exp( - PMF(r) / kT ). So exp_free_energy is g(r). This makes sense because the PMF is the energy as a function of COM distance, thus if you put the PMF in a boltzmann factor you'll get a weighting function based on the energy as a function of COM displacement, ie g(r).

        return free_energy, exp_free_energy;



    def makePlotFiles(self, bin_centers, free_energy, exp_free_energy):
        ######## PLOT SHIFTED PMF
        test_out = np.column_stack((bin_centers,free_energy))
        np.savetxt('free_energy.out',test_out)## free_energy values are now volume corrected and shfited to zero at long R. They have been shifted by the delta_PMF which is itself volume corrected.

        ######## PLOT g(r)
        test_out = np.column_stack((bin_centers,exp_free_energy))
        np.savetxt('exp_free_energy.out',test_out)



    # calculate integral of a bounded region using Simpson's Rule.
    ## Simpson's Rule requires an even number of subdivisions! (crit_bin_index must be even!)
    ## Since the index goes from 0 --> n, EVEN values of n indicate an EVEN number of subdivisions.
    def integrateSimpsons(self, r1, r2, bin_centers, function):
        Area = 0.

        r1_index = bisect_left(bin_centers,r1)## returns the index of the 'Coordinate' point in bin_centers that is closest to 'r_critical'
        r2_index = bisect_left(bin_centers,r2)## returns the index of the 'Coordinate' point in bin_centers that is closest to 'r_critical'
        index_range = r2_index - r1_index ## total number of elements in the monomer range. This will need to be an odd number for Simpson's rule.
        print 'r1 and r1_index: ',r1, r1_index
        print 'r2 and r2_index: ',r2, r2_index
        print 'bin_centers[r1_index]: ',bin_centers[r1_index]
        print 'bin_centers[r2_index]: ',bin_centers[r2_index]
        #print 'index_range = ',index_range

        #crit_bin_index_mon = index_range - 1 ## last index value of monomer_range ie. [0,monomer_range) == [0,crit_bin_index_mon]. This will need to be an even number for Simpson's rule.

        range_mod = index_range % 2 ## == 0 if divisible by 2. =/= 0 otherwise. For Simpson's Rule
        if range_mod:## This statement will execute if monomer_mod =/= 0
            index_range -= 1 ## minus instead of plus because our array doesnt have any more points in it to add.
        #    print 'Shifted the last index by 1 for Simpson\'s Rule: %s' %(crit_bin_index_mon)
        #else:## This statement will execute if monomer_mod == 0
        #    print 'Already even # of monomer subdivisions for Simpson\'s Rule: %s' %(crit_bin_index_mon)

        count = 0
        for i in range(r1_index, r2_index+1):
            #print count, i
            if count == 0 or count == (r2_index - r1_index):
                Area += function[i]
                #print 'first or last: ',i ,function[i]
            elif count%2 == 0:
                Area += 2*function[i] * bin_centers[i]**2
                #print 'even: ',i, function[i]
            else:
                Area += 4*function[i] * bin_centers[i]**2
                #print 'odd: ',i ,function[i]
            count += 1

        Area *= delta_x/3.
        
        return Area;



    def calculateKa(self, wham_file, r1, r2, max_dist):

        # read in WHAM file and assign x-axis, FE array, and dx.
        bin_centers, pmf, dx = eq.readWHAM(wham_file)

        # volume correct PMF and calculate g(r) from pmf. Add points to give value out to max_dist.
        pmf, gr = volumeCorrect_addDistance(bin_centers, pmf, dx, max_dist)

        # make files to plot PMF and g(r).
        makePlotFiles(bin_centers, pmf, gr)

        # Calculate Ka using De Jong (2011) eqn. 44

        #################################################################################
        ## NUMERATOR in De Jong eqn 44. ("dimer" radii) 

        numerator_int = integrateSimpsons(0., r1, bin_centers, gr)
        print 'Numerator Integral: ',numerator_int

        #################################################################################
        ## DENOMINATOR in De Jong eqn 44. ("monomer" radii) 

        denominator_int = integrateSimpsons(r1, r2, bin_centers, gr)
        print 'Denominator Integral: ',denominator_int

        #coeff = ( (bin_centers[-1])**3 )/((bin_centers[crit_bin_index])**3)
        coeff = (4*np.pi/3.) * ( bin_centers[-1]**3 )/( 1660. )
        print "Coefficient: ",coeff
        dA = - kT * np.log( coeff * (numerator_int / denominator_int) )
        Ka = np.exp(-dA / kT)

        #########################################################

        print "Ka = %e" %(Ka)
        Ga = -kT*np.log(Ka)
        print "Ga = %s" %(Ga)



###############################################
################ Main Program #################
###############################################
global kT
kT = 0.5921868
# max dist for 5mM is 87.2516
max_dist = 87.2516

wham_file = sys.argv[1]
r1 = float(sys.argv[2])
r2 = float(sys.argv[3])
if r2 == -1:
    r2 = max_dist


eq = Equil()

eq.calculateKa(wham_file, r1, r2, max_dist)
