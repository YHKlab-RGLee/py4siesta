#!/usr/bin/env python
#------------------------------------------------------------------------
# This script plots the E-k band diagram reading the SIESTA output file #
# systemlabel.bands. If there are bands for both spin states, band for  #
# each one are plotted                                                  #
#                                                                       #
#  Script by Javier Junquera and Pablo Aguado-Puente (2007)             #
#  Updated by Hyo Seok Kim (2013)                                       #                
#------------------------------------------------------------------------

import sys,time
#import matplotlib.pyplot as plt
import numpy as np
#from pylab import *

version = 20130115

#nt = time.localtime()
#now_time = "%s_%s_%s_%s%s%s" % (nt[0],nt[1],nt[2],nt[3],nt[4],nt[5])
#usage = "Usage: %s [systemlabel.band] [E min] [E max]" % sys.argv[0]
#foottext = '\n Thank you\n## Kim,Hyo Seok (KAIST) <softmax1986@kaist.ac.kr>'

#print "## Plotting band structure..."
#print "## Version : %s \n" % version

#Check input file

#if len(sys.argv) != 4:
#    print usage
#    print foottext
#    sys.exit(1)

#if len(sys.argv) == 4:
#    filename = str(sys.argv[1])
#    ymn = float(sys.argv[2])
#    ymx = float(sys.argv[3])


def write_oneD_one_line(file_name, array_x, array_y, label_x, label_y):
    
    oneD = open(file_name, 'w')
    if len(array_x) != len(array_y):
        raise ValueError, "array_x and array_y should have same length."
    
    head1 = "#NumField:1\n"
    head2 = "#LabelX:%s,LabelY:%s\n" % (label_x, label_y)
    oneD.write(head1); oneD.write(head2)
    head3 = "#Field1:%s,NumPoint:%i\n" % (file_name, len(array_x))
    oneD.write(head3)
    
    i_y = 0
    for x in array_x:
        line = "%12.8f    %12.8f\n" % (x, array_y[i_y]); i_y +=1
        oneD.write(line)

    oneD.close()


#Read the input file where the band structure is stored.

def DataTnBand(filename,ymn=None,ymx=None):

    f = open(filename,"r")

#=========================================Read data==============================================#

#Data_of_Band : [[Fermi level], [K-point Min, K-point Max], [Energy min, Energy Max],[The number of obitals, The number of spin components, The number of K-points]]
#xk_points : The list of all of the  k-points
#bands : the list of eigenvalnue energy
#Special_k_point : [['Gamma, "'K'"], ......]
    Data_of_Band = []
    xk_points = []
    bands = []
    special_k_point = []

#Read information of part of Data of bands. 
    for i in range(0,4):
        line = f.readline()
        list_line=line.split()
        Data_of_Band.append(list_line)

#norb is the number of orbital in the unit cell,
#nspin is the nummer of spin components,
#nkpoint is the number of k-point,
#nlines is the number of lines per kpoints containing the obital
    norb = int(Data_of_Band[3][0]) * int(Data_of_Band[3][1])
#    print "Number of orbitals in unit cell           =", norb
    nspin = int(Data_of_Band[3][1])
#    print "Number of different spin polarization     =", nspin
    nkpoint = int(Data_of_Band[3][2])
#    print "Number of k-points to represent the band  =", nkpoint
    nlines = norb / 10
    if norb < 10: nlines = 1
    nrest = norb - (nlines * 10)

#Initialize the vector hoding the bands
    for initi_band in range(0, norb):
        bands.append([])

#Read the bands and xk_points
    for initi_k in range(1, nkpoint + 1):
        for initi_lines in range(0, nlines):
            line = f.readline()
            t = line.split()
            if initi_lines == 0:
                xk_points.append(float(t[0]))
                for it in range(1, len(t)):
                    initi_band = it - 1
                    bands[initi_band].append(float(t[it])-float(Data_of_Band[0][0]))
            else:
                for it in range(0, len(t)):
                    initi_band = initi_lines * 10 + it
                    bands[initi_band].append(float(t[it])-float(Data_of_Band[0][0]))
        if nrest > 1:
            line = f.readline()
            t = line.split()
            for it in range(0,len(t)):
                initi_band = nlines*10+it
                bands[initi_band].append(float(t[it])-float(Data_of_Band[0][0]))

#Initialize the k point label
    line = f.readline()
    t = line.split()
    nlabels = int(t[0])

#Read special k_point information
    for s in range(0,int(t[0])):
        line = f.readline()
        t = line.split()
        special_k_point.append(t)


    
#====================== Plotting E-k diagram (Using Matplotlib)===========================#



#Set range of x, y axis
    Fermi_Energy = float(Data_of_Band[0][0])
    ymin = min(bands[0])
    ymax = 0
    if ymn == None and ymx == None:
	ymin = ymin-abs(Fermi_Energy-ymin) * 0.1
    	ymax = Fermi_Energy + abs(Fermi_Energy-ymin)
    else:
    	ymin = ymn
    	ymax = ymx
    xmin = float(Data_of_Band[1][0])
    xmax = float(Data_of_Band[1][1])

#Drawing structure of figure
    #fig = plt.figure()
    #fig1 = fig.add_subplot(111)
    #fig1.set_ylabel(r'$E$ $-$ $E_f$   $[eV]$', fontsize=20)
    #fig1.set_xlabel('')
    #fig1.set_xlim(xmin, xmax)
    #fig1.set_ylim(ymin, ymax)
    #xticks([])

#Plot Special K points 
    num_k = len(special_k_point)
    k_point = []
    name_k_point = []
    for i_kp in range(0, num_k):
        k_point.append(float(special_k_point[i_kp][0]))
        name_k_point.append(special_k_point[i_kp][1])
    for j_kp in range(0,num_k):
        #fig1.axvline(x=k_point[j_kp], ymin=ymin, ymax=ymax, color ='r')
        #text(k_point[j_kp], ymin-0.7, name_k_point[j_kp], size = 13, 
	#     horizontalalignment='center',verticalalignment='center')
	print j_kp+1
        write_oneD_one_line('bandSpecialKPT%3.3i.oneD' % (j_kp+1),
                            k_point[j_kp]*np.ones(100),
                            np.linspace(ymin, ymax, 100),
                            'k-vector', 'eigenvalue[eV]')
#plotting band structure(Spin consideration):plotting different color
    if nspin > 1:
        for i_plot in range(0, norb/nspin-1):
            list_X = xk_points
            list_Y = bands[i_plot]
            array_X = np.array(list_X); array_Y = np.array(list_Y)
            #fig1.plot(array_X, array_Y, color = 'black', linewidth = '4')
            write_oneD_one_line('band_%3.3i_spin1.oneD' % (i_plot+1),
                                array_X, array_Y,
                                'k-vector', 'eigenvalue[eV]')
            del list_X; del list_Y

        for i_plot in range(norb/2,norb):
            list_x = xk_points
            list_y = bands[i_plot]
            array_x = np.array(list_x); array_y = np.array(list_y)
            #fig1.plot(array_x, array_y, color = 'blue', linewidth = '1')
            write_oneD_one_line('band_%3.3i_spin2.oneD' % (i_plot-norb/2+1),
                                array_x, array_y,
                                'k-vector', 'eigenvalue[eV]')
            del list_y; del list_x
#plotting band structure(non spin consideration)
    else:
        for i_plot in range(0,norb):
            list_x = xk_points
            list_y = bands[i_plot]
            array_x = np.array(list_x); array_y = np.array(list_y)
            #fig1.plot(array_x, array_y, color = 'black', linewidth = '2')
            write_oneD_one_line('band%3.3i.oneD' % (i_plot+1),
                                array_x, array_y,
                                'k-vector', 'eigenvalue[eV]')
            del list_y; del list_x
    #savefig('%s.png' % filename)
    #plt.show()

if __name__ == '__main__':
    DataTnBand(sys.argv[1])

