#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thursday January 09 10:53:58 2020
To obtain the Electron (Charge) Density Difference
Copyright: 8th, January, Yaqiong Su, Xiamen
Y-Q.Su@tue.nl; yqsu1989@gmail.com; 398168203@qq.com
@author: Yaqiong Su, Xiamen
"""
# Ref: http://bbs.keinsci.com/thread-4476-1-1.html

import numpy as np
import datetime
import time
import sys
import linecache as lce
import pandas as pd
from scipy.integrate import simps
from operator import mul
import os

np.set_printoptions(suppress=True)
#CHGCAR0, CHGCAR1, CHGCAR2 = sys.argv[1],sys.argv[2],sys.argv[3]   # pass parameters chgcar files

######   timing   ######
start = time.time()
print '********** The averaged Chrarge Density Difference along z-axis from Yaqiong Su Xiamen **********'
print 'is getting the averaged Chrarge Density Difference along z-axis'
### current time ###
start_time = datetime.datetime.now()
print "Start time:       " + start_time.strftime('%Y.%m.%d-%H:%M:%S')   #strftime可以自定义时间的输出格式

### getting parameters such as Volume ###
CHGCAR2 = 'CHGCAR' + str(2)
locals()['CHGCAR' + str(2)] = CHGCAR2
line3 = lce.getline(CHGCAR2,3)
line4 = lce.getline(CHGCAR2,4)
line5 = lce.getline(CHGCAR2,5)
pa = float(line3.split()[0])   # obtain lattice a unit angstrom
pb = float(line4.split()[1])   # obtain lattice b unit angstrom
pc = float(line5.split()[2])   # obtain lattice c unit angstrom
Volume = pa*pb*pc       # obtain slab valumeunit A3
print 'Slab Volume is: ', Volume, ' A3'

def readzchg(i):
# read FFT-mesh of charge density from CHGCAR
    CHGCARi = 'CHGCAR' + str(i)
    locals()['CHGCAR' + str(i)] = CHGCARi
    line7 = lce.getline(CHGCARi,7)
    atoms_number = map(int,line7.split())   # map can convert string to other data-type, such as integer
    print 'atoms number:', atoms_number
#ntot_atoms = sum(int(i) for i in atoms_number)
    ntot_atoms = sum(atoms_number)
    print 'total number of atoms:', ntot_atoms
    line_FFT = ntot_atoms+10   # the line containing FFT-mesh points
    lineFFT = lce.getline(CHGCARi,line_FFT) 
    mesh = map(int,lineFFT.split()) 
    print 'FFT-mesh points:', mesh  # the mesh points along x, y, and z axis
    line_start = ntot_atoms+11 # starting line of charge density 
    xy = reduce(mul,mesh)
    inte = xy / 5  # integer
    rem = xy % 5   # reminder
    print 'integer', inte
    print 'reminder', rem
    if rem < 5 and rem > 0:
        line_end = inte + 1  + line_start - 1  # ending line of charge density 
    else:
        line_end = inte + line_start - 1  # ending line of charge density  
    chgi = 'chg' + str(i)
    locals()['chg' + str(i)] = chgi
    f1 = open(chgi,'wb')
    with open(CHGCARi,'rb') as f:
        j = 0
        while True:
            j+=1
            line = f.readline()
            if j > line_start-1 and j < line_end+1:   # the range of charge density (FFT-mesh matrix)
                f1.write(line)
            if j > line_end: 
                break
    f1.close()   # save to chg1 file

#    matrix1 = np.loadtxt(chgi)   # dtype = float (default), also can designate which columns are used
#    matrix11 = np.reshape(matrix1,(mesh[0]*mesh[1],mesh[2]))
    
    with open(chgi) as ff:
        line_N = []
        k = 0
        while True:
            line = ff.readline()
            k+=1
            if k < line_end -line_start + 1 + 1:
                line = line.strip()
                num = line.split()
                for count in range(len(num)):    
                    line_N.append(num[count])   # write element one by one
            if k > line_end + 100:
                break
    print 'total number of FFT-mesh points:', len(line_N)
    line_NN = map(float,line_N)   # convert string to float
    matrix1 = np.mat(line_NN)          
    matrix10 = np.reshape(matrix1,(mesh[2],mesh[0]*mesh[1]))
    matrix11 = matrix10.T
# The palne-averaged chrge density along z-axis 
#print matrix11.sum(axis=0) # total values of each column
    l11r = matrix11.sum(axis=0)/(mesh[0]*mesh[1]) # averged values of each column, saved to a row-list
    m11r = np.mat(l11r)  # saved to a row-matrix
    m11ci = 'm11c' + str(i)
    locals()['m11c' + str(i)] = m11ci
    m11ci = m11r.reshape(-1,1) # saved to a column-matrix
    fm = 'm11c' + str(i) + '.dat'
    locals()['m11c' + str(i) + '.dat'] = fm
    np.savetxt(fm,m11ci)
    step = 1/float(mesh[2])  # step along z-axis
    global step
    zpoints = np.arange(0,1,step) # range only constructs integer list, while np.arrange can do float
    zpc = np.mat(zpoints).reshape(-1,1) # saved to a column-matrix
    zchgi = 'zchg' + str(i)
    zchgi = np.hstack((zpc,m11ci)) # map each charge density to each z point by combining two matrix
    locals()['zchg' + str(i)] = zchgi
    file = 'zchg' + str(i) + '.dat'
    locals()['zchg' + str(i) + '.dat'] = file
    np.savetxt(file,zchgi)
    return

#readzchg(i for i in range(3))
readzchg(0)
readzchg(1)
readzchg(2)
m11c0 = np.loadtxt('m11c0.dat')
m11c1 = np.loadtxt('m11c1.dat')
m11c2 = np.loadtxt('m11c2.dat')
m11rdiff = m11c2 - m11c1 - m11c0 # get difference
m11cdiff = m11rdiff.reshape(-1,1) # row ro column   charge density multiplied Volume
m11cdiffV = m11cdiff / Volume
#step = 1/float(mesh[2])  # step along z-axis
print step
zpoints = np.arange(0,1,step) # range only constructs integer list, while np.arrange can do float
zpc = np.mat(zpoints).reshape(-1,1) # saved to a column-matrix
zchgdiff = np.hstack((zpc,m11cdiffV)) # charge density
np.savetxt('zchgdiff.dat',zchgdiff)   

os.system('rm chg0  chg1  chg2 m11c0.dat m11c1.dat  m11c2.dat  zchg0.dat  zchg1.dat  zchg2.dat')
##########   timing   #############
stop=time.time()
print("running time:     " + str(stop-start) + " seconds")
terminal_time = datetime.datetime.now()
print "Terminal time:    " + terminal_time.strftime('%Y.%m.%d-%H:%M:%S')   #strftime可以自定义时间的输出格式
print 'The averaged Chrarge Density Difference along z-axis has obtained'
print '********** The averaged Chrarge Density Difference along z-axis from Yaqiong Su Xiamen **********'
