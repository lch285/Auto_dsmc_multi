#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 13:55:58 2022

@author: luischacon
"""

from automated_multy import automated
import numpy as np
from pyDOE import lhs
import os

if __name__ == '__main__':

    """
    Basic settings
    
    """
    cluster = 'NASA' # either 'LCC' , 'NASA', or 'MCC'
    typestl = 'XRCT' # either 'XRCT' or 'Fibergen'
    
    stlfile = '200cube.stl'
    convertionfactor = '10**-6'
    
    poolsize = 6
    maxjobperrun = 25
    
    """
    Convergence test
    
    """
    convergence_flag = False
    
    if convergence_flag:
        # Just change one at a time
        
        Tconverg = 600 # C
        Pconverg = 2000 # Pa
        TandP = [Tconverg, Pconverg]
        
        if typestl == 'Fibergen':
            Length = 100 # Length of the Fibergen volume
            TandP.append(Length)
        
        spacefraction = [2] # resolution in space
        timefraction = [5, 7, 12, 15] # resolution in time
        numtimesteps = None # either a [number] or  None and fix targettimefraction and targettime accordingly
        
        
        targettimefraction = 3
        targettime = 37500
        if numtimesteps == None:
            totaltime = targettime / targettimefraction
            
            numtimesteps = np.zeros(len(timefraction))
            for num ,i in enumerate(timefraction):
                numtimesteps[num] = i*totaltime
            
        casesConverg = []
        for i in spacefraction:
            if len(timefraction) != 1:
                for num, j in enumerate(timefraction):
                    casesConverg.append([i,j,numtimesteps[num]])
            else:
                for j in timefraction:
                    for k in numtimesteps:
                        casesConverg.append([i,j,k])
    
    """ 
    Regular run
    
    """
    
    regular_flag = False
    if regular_flag:
        casesConverg = [2, 5, 62500] # Convergence criteria
        
        Temperature = np.arange( 300, 400, 100) # (start, end, incremet) end will not get executed
        Pressure = np.arange( 250, 350, 100) # (start, end, incremet) end will not get executed
        
        if typestl == 'Fibergen':
            Length = np.arange(50, 150, 100) # (start, end, incremet) end will not get executed
            
            dimension = len(Temperature)*len(Pressure)*len(Length)
            TandP = np.zeros((dimension,3))
            
            num = 0
            for T in Temperature:
                for P in Pressure:
                    for L in Length:
                        TandP[num] = [T, P, L]
                        num += 1
                        
        elif typestl == 'XRCT':
            dimension = len(Temperature)*len(Pressure)
            TandP = np.zeros((dimension,2))
            
            num = 0
            for T in Temperature:
                for P in Pressure:
                    TandP[num] = [T, P]
                    num += 1
    
    """"
    LHS run
    
    """
    
    lhs_flag = True
    
    if lhs_flag:
        casesConverg = [2, 5, 62500] # Convergence criteria
        
        Npoints = 2 # Number of points to generate
        
        Tmin = 300
        Tmax = 400
        Pmin = 1000
        Pmax = 2000
        
        if typestl == 'Fibergen':
            
            Lmin = 300
            Lmax = 600
            
            input_mtx = np.array([[Tmin, Tmax] , [Pmin, Pmax], [Lmin, Lmax]])
            dim = 3
            #Variable Matrix
            TandP=np.zeros((Npoints,dim))
            LHD = lhs(dim,samples=Npoints,criterion="m")
            for i in range(Npoints):
                for j in range(dim):
                    TandP[i,j] = input_mtx[j,0] + LHD[i,j]*(input_mtx[j,1] - input_mtx[j,0])
                    
        elif typestl == 'XRCT':
            input_mtx = np.array([[Tmin, Tmax] , [Pmin, Pmax]])
            dim = 2
            #Variable Matrix
            TandP=np.zeros((Npoints,dim))
            LHD = lhs(dim,samples=Npoints,criterion="m")
            for i in range(Npoints):
                for j in range(dim):
                    TandP[i,j] = input_mtx[j,0] + LHD[i,j]*(input_mtx[j,1] - input_mtx[j,0])
            
    
    """
    Start of the automated script
    """
    
    
    if (int(lhs_flag) + int(regular_flag) + int(convergence_flag)) == 1:
        
        if os.path.exists(stlfile) or typestl == 'Fibergen':
        
            automated(cluster, typestl, stlfile, convertionfactor, convergence_flag , TandP, casesConverg, poolsize, maxjobperrun)
            
        else:
            print('stl file does not exist!!')
    
    
    
    

