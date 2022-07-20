# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 15:54:51 2022

@author: luisa
"""

import time
import numpy as np
import multiprocessing as mp
import math
import os


def process(queue, results):   #(setlines,micro_domain,variables,di):
    
    setlines,micro_domain,variables,di = queue.get()
    variable_mtx=np.zeros((len(setlines),40))
    ComputedValues = np.zeros((9))
    
    # input("Press Enter to continue...")
    for no, line in enumerate(setlines):
        s = tuple(line.split())
        n=len(s)
        for i in range(1,n):
            variable_mtx[no,i-1]=float(s[i])
            
    # print(variable_mtx)
    variable_mtx=variable_mtx.transpose()
    variable_mtx = variable_mtx[~np.all(variable_mtx == 0, axis=1)]
    variable_mtx=variable_mtx.transpose()
    # print(variable_mtx)
    
    ComputedValues= ComputeVal(variable_mtx, micro_domain, variables, di)
    
    results.put( ComputedValues)
    # return ComputedValues
    # print(ComputedValues)
    
    # input("Press Enter to continue...")

def ComputeVal(variable_mtx,micro_domain,variables,di):

    T1=0
    T2=0
    Rho1=0
    Rho2=0
    u1=0
    u2=0
    left=0
    right=0
    integral=0
    
    for i in range (len(variable_mtx)):
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
            
    return T1,Rho1,u1,left,T2,Rho2,u2,right,integral



def pp_parallel_fast(temp_number,domain_extend,pathmain):
    os.chdir(pathmain+'/Results_multi/dsmc_temp%d' %temp_number)
    start = time.time()
    
    nproces = mp.cpu_count()
    
    
    
    
    chunksize = 2*10**6
    
    filename = 'flow.output'
    
    variables=tuple()
    species=tuple()
    micro_domain=np.zeros((1,2))
    f1=open('dsmc.input','r')
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
            

    micro_domain[0,0]=sep_lim[0,0]+domain_extend # if not a voxell file uncoment
    micro_domain[0,1]=sep_lim[0,1]-domain_extend # if not a voxell file uncoment
    f1.close()
    
    #Permeability Force
    f=open('collision.list','r')
    for line in f:
        if species[0] in line:
            s=tuple(line.split())
            dref=float(s[1])
            omega=float(s[2])
            Tref=float(s[3])
            break
    f.close()
    f=open('species.list','r')
    for line in f:
        if species[0] in line:
            s=tuple(line.split())
            Molmass=float(s[2])
            break
    f.close()
    
    
    end=time.time()
    print('Get properties: %.4f' % (end-start))
    
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
    
    
    with open(filename, mode="r", encoding="utf8") as file_obj:
        
        lines=[]
        proceses =[]
        queue = mp.Queue()
        results = mp.Queue()
        
        
        pool = mp.Pool(nproces)
        
        end0 = time.time()
        end1 = time.time()
        
        j = 1
        flag = 0
        flag_ncell = 0 
        for num, line in enumerate(file_obj):
            
            
            if line == timestep:
                end = time.time()
                print('Found first line: %.4f' % (end-end0))
                
                flag = 1
                firstline = num
             
                
            if flag_ncell:
                act_ncell = int(line)
                flag_ncell =0
                
            if 'ITEM: NUMBER OF CELLS' in line:
                flag_ncell = 1

            if flag and 'ITEM: TIMESTEP' in line:
                print('loop broken!!!')
                num=num-1
                break
            if (flag and num >= firstline+8):
    
                lines.append(line)
                
                if num-firstline-8 >= chunksize*j:
                    end = time.time()
                    
                    percentCompleted = 100 - 100* np.divide( np.abs( np.subtract( chunksize*j,act_ncell)),act_ncell) 
                    
                    print('Entering extraction Group %i (%.2f' % (j,percentCompleted) +' %):')
                    # results.append( pool.starmap( process,[(lines, micro_domain, variables, di)]))
                    queue.put((lines, micro_domain, variables, di))
                    p = mp.Process(target=(process), args=(queue,results))
                    proceses.append(p)
                    p.start()
                    
                    
                    end0 = time.time()
                    print('Exiting extraction Group %i after: %.4f'% (j, end0-end) )
                    
                    j += 1
                    lines = []
                elif num-firstline-8 +1 == act_ncell:
                    end = time.time()
                    
                    percentCompleted = 100
                    print('Entering extraction Group %i (%.2f' % (j,percentCompleted) +' %):')
                    # results.append( pool.starmap( process,[(lines, micro_domain, variables, di)]))
                    queue.put((lines, micro_domain, variables, di))
                    p = mp.Process(target=(process), args=(queue,results))
                    proceses.append(p)
                    p.start()
                    
                    
                    end0 = time.time()
                    print('Exiting extraction Group %i after: %.4f'% (j, end0-end) )
                    
                    
                
    print('Waiting for all resutls to compute!')      
    for p in proceses:
        p.join()
    
    end = time.time()
    print('Finish Computation: %.4f' % (end - end1) )
    
    ComputedValues = np.zeros((j,9))
    for i, p in enumerate(proceses):
        if i > j:
            break
        ComputedValues[i,:] = results.get()
            
    SumedValues = np.zeros(9)
    for i in range(j):
        SumedValues += ComputedValues[i]

    T1,Rho1,u1,left,T2,Rho2,u2,right,integral = SumedValues


    T1=np.divide(T1,left)
    Rho1=np.divide(Rho1,left)
    u1=np.divide(u1,left)
    T2=np.divide(T2,right)
    Rho2=np.divide(Rho2,right)
    u2=np.divide(u2,right)
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
        
    end0=time.time()
    print('Finish Calculations: %.4f' % (end0-end))
    
    # print(viscosity)
    # print(P1)
    # print(P2)
    # print((P1+P2)/2)
    # print(u1)
    # print(u2)
    # print(Mass_rate1)
    # print(Mass_Rate2) 
    # print(Perm_force)
    # print(T)
    # print(T1)
    # print(T2)
    # print(K)
    # print('F=',Perm_force)
    print('00 Ar 000',domain_extend*10**6, T, (P1+P2)/2, K, Perm_force)
    # print('Chunks size:', chunksize)
    # print('With', nproces, 'processes.')
        
        
    end = time.time()
    
    print('Total time: %.4f' % (end-start))
    return T, (P1+P2)/2, K , Perm_force
