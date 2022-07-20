#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 20:33:03 2021

@author: vijaybmohan
"""

import os
from evtk.hl import structuredToVTK
import numpy as np
import math
import csv
def postprocess(temp_number,domain_extend,pathmain):
    os.chdir(pathmain+'/dsmc_temp%d' %temp_number)
    timestep='10000\n'
    time_flag=0
    f=open(pathmain+'/dsmc_temp%d/flow.output' %temp_number,'r')
    for num, line in enumerate(f, 1):
        if line==timestep:
            first_line=num
            time_flag=1
        if time_flag:
            if 'ITEM: TIMESTEP' in line:
                num=num-1
                break
    print(num)
    print(first_line)
    f.close()
    variables=tuple()
    species=tuple()
    micro_domain=np.zeros((1,2))
    f1=open(pathmain+'/dsmc_temp%d/dsmc.input' %temp_number,'r')
    for line in f1:
        if ('dimension') in line:
            di=int(line[10])
            cells_di=np.zeros((di,))
            sep_di=np.zeros((di,))
            sep_lim=np.zeros((di,2))
        if ('create_box') in line:
            s=tuple(line.split())
            for i in range(di):
                sep_di[i]=float(s[2*i+2])-float(s[2*(i-1)+3])
                sep_lim[i,0]=float(s[2*(i-1)+3])
                sep_lim[i,1]=float(s[2*i+2])
        if ('create_grid') in line:
            s=tuple(line.split())
            for i in range(di):
                cells_di[i]=s[i+1]    
        if ('region region_microstructure block') in line:
            s=tuple(line.split())
            micro_domain[0,0]=float(s[3])
            micro_domain[0,1]=float(s[4])
        if ('species species.list') in line:
            s=tuple(line.split())
            n=s.index('species.list')
            species=species+s[n+1:]
        if ('compute') in line:
            if 'species' in line:
                s=tuple(line.split())
                n=s.index('species')
                variables=variables+s[n+1:]
        if ('run') in line:
            s=tuple(line.split())
            timestep = s[1]+'\n'
            timefloat = float(s[1])
    micro_domain[0,0]=sep_lim[0,0]+domain_extend
    micro_domain[0,1]=sep_lim[0,1]-domain_extend
    f1.close()
        
    # Dimensions
    nx, ny = int(cells_di[0]), int(cells_di[1])
    lx, ly = sep_di[0], sep_di[1]
    dx, dy = lx/nx, ly/ny
    ncells = nx * ny
    if di==3:
        nz=int(cells_di[2])
        lz=sep_di[2]
        dz=lz/nz
        ncells=ncells*nz
        
    # Variables
    variable_mtx=np.zeros((num-first_line-7,20))
    f=open(pathmain+'/dsmc_temp%d/flow.output' %temp_number,'r')
    for no,line in enumerate(f,1):
        if no>=first_line+8:
            #print(no)
            s=tuple(line.split())
            n=len(s)
            for i in range(1,n):
                variable_mtx[no-first_line-8,i-1]=float(s[i])
        if no==num:
            break
            
    f.close()
        
    #Permeability Force
    f=open(pathmain+'/dsmc_temp%d/collision.list' %temp_number,'r')
    for line in f:
        if species[0] in line:
            s=tuple(line.split())
            dref=float(s[1])
            omega=float(s[2])
            Tref=float(s[3])
            break
    f.close()
    f=open(pathmain+'/dsmc_temp%d/species.list' %temp_number,'r')
    for line in f:
        if species[0] in line:
            s=tuple(line.split())
            Molmass=float(s[2])
            break
    f.close()
    variable_mtx=variable_mtx.transpose()
    variable_mtx = variable_mtx[~np.all(variable_mtx == 0, axis=1)]
    variable_mtx=variable_mtx.transpose()
    T1=0
    T2=0
    Rho1=0
    Rho2=0
    u1=0
    u2=0
    left=0
    right=0
    integral=0
    for i in range(num-first_line-7):
        if variable_mtx[i,0]<micro_domain[0,0]:
            n=variables.index('temp')
            T1=T1+variable_mtx[i,2+2*di+n-1]
            n=variables.index('massrho')
            Rho1=Rho1+variable_mtx[i,2+2*di+n-1]
            n=variables.index('u')
            u1=u1+variable_mtx[i,2+2*di+n-1]
            left+=1
        if variable_mtx[i,3]>micro_domain[0,1]:
            n=variables.index('temp')
            T2=T2+variable_mtx[i,2+2*di+n-1]
            n=variables.index('massrho')
            Rho2=Rho2+variable_mtx[i,2+2*di+n-1]
            n=variables.index('u')
            u2=u2+variable_mtx[i,2+2*di+n-1]
            right+=1
        if variable_mtx[i,0]>micro_domain[0,0] and variable_mtx[i,3]<micro_domain[0,1]:
            n=variables.index('u')
            integral=integral+(variable_mtx[i,6]*variable_mtx[i,2+2*di+n-1])
    T1=T1/left
    Rho1=Rho1/left
    u1=u1/left
    T2=T2/right
    Rho2=Rho2/right
    u2=u2/right
    T=(T1+T2)/2
    P1=Rho1*8.314*T1/(Molmass*6.022*10**23)
    P2=Rho2*8.314*T2/(Molmass*6.022*10**23)
    if di==2:
        Mass_rate1=Rho1*u1*ly
        Mass_Rate2=Rho2*u2*ly
    if di==3:
        Mass_rate1=Rho1*u1*ly*lz
        Mass_Rate2=Rho2*u2*ly*lz
    Mass_rate=(Mass_rate1+Mass_Rate2)/2
    viscosity=((15*math.sqrt(math.pi*Molmass*1.380649*10**-23*Tref))/(2*(5-2*omega)*(7-2*omega)*math.pi*dref**2))*((T/Tref)**omega)
    l_sample=abs(micro_domain[0,1]-micro_domain[0,0])
    if di==2:
        Perm_force=(viscosity*Mass_rate*T*8.314*l_sample)/(ly*Molmass*6.022*10**23*(P1-P2))
        K=(integral*viscosity)/(ly*(P1-P2))
    if di==3:
        Perm_force=(viscosity*Mass_rate*T*8.314*l_sample)/(ly*lz*Molmass*6.022*10**23*(P1-P2))
        K=(integral*viscosity)/(ly*lz*(P1-P2))
    return T,(P1+P2)/2,K, Perm_force,  timefloat  
            
