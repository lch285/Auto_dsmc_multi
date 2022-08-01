#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 12:51:14 2021

@author: vijaybmohan
"""

import numpy as np
import os
import shutil
import fileinput
import sys
import time
import math
import subprocess
from pp_generalized_auto import postprocess


def loop_process(cluster, typestl, stlfile, convertionfactor, x,caseConverg, convergence_flag, ident,sims,pathmain,MainName):
    temp_number=ident[0,0]+sims
    #Create temp folder
    path=pathmain+'/dsmc_temp%d' %temp_number
    
    try:
       if not os.path.exists(path):
           os.mkdir(path)
    except OSError as err:
       print(err) 
    source_dir = pathmain+MainName
    target_dir = pathmain+'/dsmc_temp%d' %temp_number
      
    file_names = os.listdir(source_dir)
        
    for file_name in file_names:
        if os.path.isfile(file_name):
            shutil.copy(os.path.join(source_dir, file_name), target_dir)
    
    os.chdir(pathmain+'/dsmc_temp%d' %temp_number)
    file_names = os.listdir(target_dir)       
    
    for filename in file_names:
        if filename.endswith('.out'):
            os.unlink(filename)
    
    if typestl == 'Fibergen' and not convergence_flag:
        # Model for Porosity of FiberForm
        a =      0.9042  
        b =  -6.574e-05  
        c =     -0.4612  
        d =    -0.01315  
        porosity = a*np.exp(b*x[2]) + c*np.exp(d*x[2])
        
        if x[2] > 400:
            porosity = 0.8783369779174772

        
     
        #Modify inputdeck.in
        f_in='inputdeck.in'
        for line in fileinput.input(f_in,inplace=1):
            if 'eps_nom= porosity' in line:
                line=line.replace('porosity',str(porosity),1)
            if 'x_box= length' in line:
                line=line.replace('length',str(x[2])+'d-6',1)
            if 'y_box= length' in line:
                line=line.replace('length',str(x[2])+'d-6',1)
            if 'z_box= length' in line:
                line=line.replace('length',str(x[2])+'d-6',1)
            sys.stdout.write(line)
            
     
        #Create grid.stl
        n=0
        n1=1
        os.system('/project/sjpo228_uksr/LuisChacon/git/fibergen/fibergen')
        path_fibergen=(pathmain+'/dsmc_temp%d/microstructure_values.dat' %temp_number)
        flag_fiber_stop=1
        while(flag_fiber_stop):
            if os.path.isfile(path_fibergen):
                n=os.path.getsize(path_fibergen)
                time.sleep(2)
                n1=os.path.getsize(path_fibergen)
            else:
                time.sleep(2)
                print('creating',flush=True)
            if n==n1:
                flag_fiber_stop=0
        
        #File conversion (grid.stl--fibergen.sparta)
        subprocess.Popen(["python", "stl2surf.py","grid_physical.stl","fibergen.sparta"])
        path_sparta=(pathmain+'/dsmc_temp%d/fibergen.sparta' %temp_number)
        flag_sparta=1
        n=0
        n1=1
        while(flag_sparta):
            if os.path.isfile(path_sparta):
                n=os.path.getsize(path_sparta)
                time.sleep(2)
                n1=os.path.getsize(path_sparta)
            else:
                time.sleep(2)
                print('converting',flush=True)
            if n==n1:
                flag_sparta=0
                
    elif typestl == 'XRCT':
        porosity = 00
        #File conversion (grid.stl--fibergen.sparta)
        
        path_sparta=(pathmain+'/dsmc_temp%d/fibergen.sparta' %temp_number)
        
        if os.path.exists(path_sparta):
            os.remove(path_sparta)
        
        subprocess.Popen(["python", "stl2surf.py", stlfile ,"fibergen.sparta", convertionfactor])
        
        flag_sparta=1
        n=0
        n1=1
        while(flag_sparta):
            if os.path.isfile(path_sparta):
                n=os.path.getsize(path_sparta)
                time.sleep(2)
                n1=os.path.getsize(path_sparta)
            else:
                time.sleep(2)
                print('converting',flush=True)
            if n==n1:
                flag_sparta=0
        
        path_fibergen=(pathmain+'/dsmc_temp%d/microstructure_values.dat' %temp_number)
        
        if os.path.exists(path_fibergen):
            os.remove(path_fibergen)
        
        with open(path_sparta, 'r') as f:

            firstline_flag = 0 
            lastline_flag = 0
            
            npoint = 0 
            xminstl = 0
            xmaxstl = 0
            yminstl = 0
            ymaxstl = 0
            zminstl = 0
            zmaxstl = 0
                
            for num, line in enumerate(f):
                if line == 'Points\n':
                    firstlinenum = num
                    firstline_flag =1
                
                if firstline_flag and num >= firstlinenum+2:
                    if line == '\n':
                        lastline_flag = 1
                        
                    if lastline_flag == 0:
                        s=tuple()
                        s = line.split()
                        for i in range(len(s)):
                            s[i] = float(s[i])
                            
                        if npoint == 0:

                            xminstl = s[1]
                            xmaxstl = s[1]
                            yminstl = s[2]
                            ymaxstl = s[2]
                            zminstl = s[3]
                            zmaxstl = s[3]                    
                            npoint += 1
                        
                        if s[1] < xminstl:
                            xminstl = s[1]                 
                        if s[1] > xmaxstl:
                            xmaxstl = s[1]             
                        if s[2] < yminstl:
                            yminstl = s[2]                   
                        if s[2] > ymaxstl:
                            ymaxstl = s[2]               
                        if s[3] < zminstl:
                            zminstl = s[3]                    
                        if s[3] > zmaxstl:
                            zmaxstl = s[3]
                            
                        npoint += 1

        with open(path_fibergen, 'w') as fw:
            fw.write('%s %s %s %s %s %s %s 1.0' % (porosity, xminstl, xmaxstl, yminstl, ymaxstl, zminstl, zmaxstl) )
        
            
    #Choose species
    # species_list=['CO','N2','Ar','O2','CO2']
    # random.seed(random.random())
    # spec=random.choice(species_list)
    spec='Ar'
    f='dsmc.input'
    for line in fileinput.input(f,inplace=1):
        if 'species species.list' in line:
            line=line.replace('gas_species',spec,1)
        if 'mixture inflowgas' in line:
            line=line.replace('gas_species',spec,1)
        sys.stdout.write(line)
    
    #Calculate mean free path and collision time
    f5=open(path_fibergen,'r')
    for line in f5:
        s1=tuple(line.split())
    f5.close()
    porosity=float(s1[0])
    domain_extend= (float(s1[2]) - float(s1[1]))
    xmin=float(s1[1])-domain_extend
    xmax=float(s1[2])+domain_extend
    sideextend=domain_extend*0.005
    ymin=float(s1[3])-sideextend
    ymax=float(s1[4])+sideextend
    species=tuple()
    f6=open(pathmain+'/dsmc_temp%d/dsmc.input' %temp_number,'r')
    for line in f6:
        if ('dimension') in line:
            di=int(line[10])
        if ('species species.list') in line:
            s=tuple(line.split())
            n=s.index('species.list')
            species=species+s[n+1:]
    f6.close()
    f7=open(pathmain+'/dsmc_temp%d/collision.list' %temp_number,'r')
    for line in f7:
        if species[0] in line:
            s=tuple(line.split())
            dref=float(s[1])
            omega=float(s[2])
            Tref=float(s[3])
            break
    f7.close()
    f8=open(pathmain+'/dsmc_temp%d/species.list' %temp_number,'r')
    for line in f8:
        if species[0] in line:
            s=tuple(line.split())
            Molmass=float(s[2])
            break
    f8.close()    
    particle_density=(x[1]*6.023*10**23)/(8.314*x[0])
    mean_free_path=(1/(math.sqrt(2)*particle_density*math.pi*dref**2))*((x[0]/Tref)**(omega-0.5))
    collision_time=(math.sqrt((Molmass*6.023*10**23)/(4*math.pi*8.314*Tref)))*(1/(2*particle_density*dref**2))*((Tref/x[0])**(1-omega))
    timefraction = caseConverg[1]
    time_step=collision_time/timefraction
    spacefraction = caseConverg[0]
    xcells=int((xmax-xmin)*spacefraction/mean_free_path)
    ycells=int((ymax-ymin)*spacefraction/mean_free_path)
    if xcells<201:
        xcells=201
    if ycells<67:
        ycells=67
    if di == 3:
        zmin=float(s1[5])-sideextend
        zmax=float(s1[6])+sideextend
        total_vol=(xmax-xmin)*(ymax-ymin)*(zmax-zmin)
        zcells=int((zmax-zmin)*spacefraction/mean_free_path)
        if zcells<67:
            zcells=67
        ncells=xcells*ycells*zcells
        vol_cell=total_vol/ncells
    else:
        total_vol=(xmax-xmin)*(ymax-ymin)
        ncells=xcells*ycells
        vol_cell=total_vol/ncells
    particle_count=vol_cell*particle_density
    particle_ratio=particle_count/50
    grid_cut=2.5*mean_free_path
    
    #Modify input file
    f='dsmc.input'
    for line in fileinput.input(f,inplace=1):
        if 'create_box' in line:
            line=line.replace('xmin',str(xmin),1)
            line=line.replace('xmax',str(xmax),1)
            line=line.replace('ymin',str(ymin),1)
            line=line.replace('ymax',str(ymax),1)
            line=line.replace('zmin',str(zmin),1)
            line=line.replace('zmax',str(zmax),1)
        if 'global fnum' in line:
            line=line.replace('particle_ratio',str(particle_ratio),1)
            line=line.replace('particle_density',str(particle_density),1)
            line=line.replace('temp_replace',str(x[0]),1)
            line=line.replace('grid_cut',str(grid_cut),1)
        if 'create_grid' in line:
            line=line.replace('xcells',str(xcells),1)
            line=line.replace('ycells',str(ycells),1)
            if di == 3:
                line=line.replace('zcells',str(zcells),1)
        if 'timestep' in line:
            line=line.replace('time_step',str(time_step),1)
        if 'mixture inflowgas' in line:
            line=line.replace('temp_replace',str(x[0]),1)
        if 'surf_collide 1 diffuse' in line:
            line=line.replace('temp_replace',str(x[0]),1)
        if 'fix in emit/face' in line:
            line=line.replace('P1',str(x[1]+50),1)
            line=line.replace('temp_replace',str(x[0]),1)
        if 'fix out emit/face' in line:
            line=line.replace('P2',str(x[1]-50),1)
        if 'avgtime' in line:
            line=line.replace('avgtime',str(int(caseConverg[2]/2)),1)
        if 'totaltimesteps' in line:
            line=line.replace('totaltimesteps',str(caseConverg[2]),1)
        sys.stdout.write(line)
    
    # Get number of nodes and processes


    if cluster == 'NASA':
        N_processors = ncells/50000
        N_processors_node = 20
        N_nodes = int(N_processors/N_processors_node+1)
        total_processors = N_processors_node*N_nodes
        if total_processors > 20000:
            total_processors = int(20000/N_processors_node)*N_processors_node
            N_nodes = total_processors/N_processors_node
        
        submitcommand = 'qsub'
        f3='submitNAS.sh' 
        for line in fileinput.input(f3,inplace=1):
            line=line.replace('nodes:ncpus=processors:mpiprocs=processors','%i:ncpus=%i:mpiprocs=%i'% (N_nodes, N_processors_node, N_processors_node),1)
            line=line.replace('slurmtemp','slurm%i'%temp_number,1)
            line=line.replace('totalprocessors','%i'%total_processors,1)
            sys.stdout.write(line)
        
    elif cluster == 'LCC':
        N_processors = ncells/60000
        N_processors_node = 48
        N_nodes = int(N_processors/N_processors_node+1)
        total_processors = N_processors_node*N_nodes
        if total_processors > 480:
            total_processors = int(480/N_processors_node)*N_processors_node
            N_nodes = total_processors/N_processors_node
        
        submitcommand = 'sbatch'
        f3='submitLCC.sh' 
        for line in fileinput.input(f3,inplace=1):
            line=line.replace('#SBATCH --job-name=np_101','#SBATCH --job-name=np_10%d' % temp_number,1)
            line=line.replace('#SBATCH -N nodes','#SBATCH -N %i'% N_nodes,1)
            line=line.replace('totalprocessors','%i'% total_processors,1)
            line=line.replace('totaltimesteps','%i'% caseConverg[2],1)
            sys.stdout.write(line)

    elif cluster == 'MCC':
        N_processors = ncells/70000
        N_processors_node = 128
        N_nodes = int(N_processors/N_processors_node+1)
        total_processors = N_processors_node*N_nodes
        if total_processors > 1408:
            total_processors = int(1408/N_processors_node)*N_processors_node
            N_nodes = total_processors/N_processors_node

        submitcommand = 'sbatch'
        f3='submitMCC.sh' #change to submit file
        for line in fileinput.input(f3,inplace=1):
            line=line.replace('#SBATCH --job-name=np_101','#SBATCH --job-name=np_10%d' % temp_number,1)
            line=line.replace('--nodes=nodes','--nodes=%i'%N_nodes,1)
            line=line.replace('totalprocessors','%i'%total_processors,1)
            sys.stdout.write(line)

        

    
    #Keep track of all simulations
    path_member_log=pathmain+MainName
    member_log=os.path.join(path_member_log,'member_log.txt' )
    f_member=open(member_log,'a')
    if os.path.getsize(member_log)==0:
        f_member.write('Temp Gas  spaceReso  timeReso Ntimesteps  Porosity  Length_Scale  Average_Temp  Average_Pressure  Eff_Permeability  Permeability_Force\n')
        f_member.close()
    
    os.system('%s %s' % (submitcommand, f3))
    
    #check for stop
    pathre=pathmain+'/Results_multi/dsmc_temp%d' %(temp_number)
    try:
       if not os.path.exists(pathre):
           os.mkdir(pathre)
    except OSError as err:
       print(err)

    flag_stop=1
    inf_count = 0
    foundlast_lines = 0
    ERROR_flag = 0
    identnum = 'lol'
    while flag_stop:
        file_names=os.listdir(target_dir) # pathmain+'/dsmc_temp%d' %temp_number
        for file_name in file_names:
            
            if ERROR_flag and file_name.endswith(identnum):
                
                os.remove(target_dir+file_name)
                ERROR_flag = 0
                foundlast_lines = 0
            else:
                foundlast_lines = 0 
                if cluster == 'LCC' or cluster == 'MCC':
                    if file_name.endswith('.out'):
                        if 'slurm' in file_name:
                            time.sleep(3)
                            file_check=(os.path.join(target_dir,file_name))
                            with open(file_check, 'r') as f4:
                                last_lines = f4.readlines()[-5:]
                                foundlast_lines = 1
                            break
                elif cluster == 'NASA':
                    if 'slurm' in file_name:
                        time.sleep(3)
                        file_check=(os.path.join(target_dir,file_name))
                        with open(file_check, 'r') as f4:
                            last_lines = f4.readlines()[-5:]
                            foundlast_lines = 1
                        break
                
        if foundlast_lines:
            for last_line in last_lines:
                if 'Histogram:' in last_line:
                    print('Found Histogram in temp %i' % temp_number)
                    time.sleep(3)
                    flag_stop=0
                    shutil.copy(file_check, os.path.join(pathre,file_name))
                    break
                
                elif 'ERROR' in last_line and 'UCX  ERROR' not in last_line:

                    inf_count +=1
                    print('Trying to fix ERROR in temp%d: %i' % (temp_number, inf_count))
                    if inf_count >5:
                        break
                    
                    xcellsold = xcells
                    ycellsold = ycells
                    zcellsold = zcells
                    xcells += 10 
                    ycells = int(xcells/3)
                    zcells = int(xcells/3)
                    
                    
                    particle_ratioold = particle_ratio
                    ncells =xcells*ycells*zcells
                    vol_cell=total_vol/ncells
                    particle_count=vol_cell*particle_density
                    particle_ratio=particle_count/50
                    
                    
                    #Modify input file
                    f='dsmc.input'
                    for line in fileinput.input(f,inplace=1):
                        if 'global fnum' in line:
                            line=line.replace('%s'%particle_ratioold,str(particle_ratio),1)
        
                        if 'create_grid' in line:
                            line=line.replace('%s'% xcellsold ,str(xcells),1)
                            line=line.replace('%s'% ycellsold,str(ycells),1)
                            if di == 3:
                                line=line.replace('%s'% zcellsold,str(zcells),1)
                        sys.stdout.write(line)
                    index = file_name.rfind(".")
                    identnum = file_name[index+1:]
                    os.remove(target_dir+file_name)
                    os.system('%s %s' % (submitcommand,f3))
                    ERROR_flag = 1
                elif 'UCX  ERROR' in last_line:
                    print('Job does not have enough memory!')
                    flag_stop = 0    
                    
                    
                else:
                    print('DSMC running %d' %temp_number, flush=True)
                    time.sleep(3)
                    
            if inf_count > 5:
                flag_stop = 0
                break
        
        else:
            print('DSMC running %d' %temp_number, flush=True)
            time.sleep(3)
            
    path_output=(pathmain+'/dsmc_temp%d/flow.output' %temp_number)
    flag_out_stop=1
    n=0
    n=1
    while(flag_out_stop):
        if os.path.isfile(path_output):
            n=os.path.getsize(path_output)
            time.sleep(2)
            n1=os.path.getsize(path_output)
        else:
            time.sleep(2)
            print('creating',flush=True)
        if n==n1:
            flag_out_stop=0
            
    y_out = [spacefraction,timefraction, None ]
    
    x_out=np.zeros((1,6))
    if ncells <= 30000000:
        print('Entering series post-process in temp%i'%temp_number)
        #Postprocessing
        x_out[0,0]=porosity
        x_out[0,1]=domain_extend*10**6 # lengthScale
        [x_out[0,2],x_out[0,3],x_out[0,4],x_out[0,5], y_out[2] ]=postprocess(temp_number,domain_extend,pathmain)  # returning T, (P1+P2)/2, K, Perm_force,  timefloat  
        
    else:
        print('File to big added to parallelcases, temp%i'%temp_number)
        # add to bigcases.txt the temp_number
        with open(pathmain+MainName+'/bigcases.txt','a') as f:
            f.write('%i\n' % temp_number)
        
        f1=open(pathmain+'/dsmc_temp%d/dsmc.input' %temp_number,'r')
        for line in f1:
            if ('run') in line:
                s=tuple(line.split())
                timefloat = float(s[1])
        y_out[2] = timefloat      
        x_out[0,0] =  porosity
        x_out[0,1] = domain_extend*10**6 
        x_out[0,2] =   temp_number
        
        
    
    z_out = np.hstack((y_out, x_out[0])) # 2  1  50  0.770525  100  597.692  1984  1.45948644e-09  7.90595398e-06
    #Simulations results
     
    #copy files to results folder
    files2save = ['dsmc.input', 'collision.list', 'species.list', 'inputdeck.in', 'flow.output', 'DSMC_script_multy.py', \
                 'microstructure_values.dat', 'pp_parallel_auto.py', 'pp_parallelFast.py', 'settings.py']
        
    for file in files2save:
        savepath = '/dsmc_temp%d/%s' % (temp_number, file)
        shutil.copy(pathmain+ savepath ,pathmain+'/Results_multi' + savepath)

    
    pathf=os.path.join(pathmain+'/Results_multi/dsmc_temp%d' %(temp_number),'log.txt')
    
    print('z_out is', z_out)
    if x_out[0,3] != 0:
        f_member=open(member_log,'a')
        with open (pathf,'w') as f_log:
            for j in range(0,len(z_out)-1):
                if j==0:
                    
                    f_log.write('%s ' %spec)
                    f_member.write(' %s    ' %temp_number)
                    f_member.write('%s    ' %spec)
                else:
                    f_log.write('%0.3f  ' %z_out[j-1])
                    f_member.write('%0.4f    ' %z_out[j-1])
                    
            f_log.write('%s' %z_out[-2])
            f_member.write('%s    ' %z_out[-2])
            
            f_log.write('%s' %z_out[-1])
            f_member.write('%s\n' %z_out[-1])
        f_member.close()
    
    os.chdir(pathmain+MainName)
    
    #remove temp dsmc folder
    if os.path.isdir(path):
        shutil.rmtree(path, ignore_errors=True)
    print('Returning from temp %i!'%temp_number)
    
    return z_out
