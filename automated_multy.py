#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 21:23:45 2021

@author: vijaybmohan
"""
import fileinput
import sys
import time
import os
import numpy as np
import multiprocessing as mp
from DSMC_script_multy import loop_process
from pp_parallel_auto import pp_parallel_fast
from alert import alert

def automated(cluster, typestl, stlfile, convertionfactor, convergence_flag , TandP, casesConverg, poolsize, maxjobperrun, species_var):
    
    # Get path and main folder name
    path = os.getcwd()
    print(path)
    index = path.rfind("/")
    pathmain = path[:index]
    print(pathmain)
    MainName = path[index:]
    print(MainName)
    
    
    #Create results directory
    pathr=pathmain + '/Results_multi'
    try:
        if not os.path.exists(pathr):
            os.mkdir(pathr)
    except OSError as err:
            print(err)
    
    #Read Input
    s_re=tuple()
    f_re=open('Restart.txt','r')
    for num, line in enumerate(f_re, 1):
        if num==1:
            run_no=int(line)
            if run_no==1:
                break
        else:
            if num==2:
                s_re=line.split()
                n_re=len(s_re)
                restart_mtx=np.zeros((1,n_re))
                for i in range(n_re):
                    restart_mtx[0,i]=float(s_re[i])
            else:
                s_re=line.split()
                n_re=len(s_re)
                temp_mtx=np.zeros((1,n_re))
                for i in range(n_re):
                    temp_mtx[0,i]=float(s_re[i])
                restart_mtx=np.row_stack((restart_mtx,temp_mtx))
    f_re.close()
    
    if convergence_flag:

        variable_mtx = casesConverg
        input_files = len(casesConverg)
        num_lines =  len(casesConverg[0])
        sims=0
        
    else:
        if run_no==1:  
            sims=0
            variable_mtx = TandP
            input_files = len(TandP)
            num_lines =  len(TandP[0])
            
            # f_input=open('textinput.txt','r')
            # s=tuple()
            # num_lines = sum(1 for line in open('textinput.txt'))
            # input_mtx=np.zeros((num_lines,3))
            # i=0
            # input_files=1
            # dim_mtx=np.zeros((num_lines,1))
            # for line in f_input:
            #     s=line.split()
            #     for j in range(3):
            #         input_mtx[i,j]=float(s[j])
            #     dim_mtx[i]=(input_mtx[i,2]-input_mtx[i,0])/input_mtx[i,1]
            #     input_files=int(input_files*dim_mtx[i])
            #     i+=1
            # f_input.close()
            
            # #Variable Matrix
            # variable_mtx=np.zeros((input_files,num_lines))  
            # skip=1
            # for j in range(num_lines-1,-1,-1):
            #     for i in range(input_files):
            #         variable_mtx[i,j]=input_mtx[j,0]+((int(i/skip))%dim_mtx[j])*input_mtx[j,1]
            #     skip=skip*dim_mtx[j]
            
            
        else:
            s_out=tuple()
            flag=0
            f_out=open('member_log.txt','r')
            for num, line in enumerate(f_out, 1):
                if num>1:
                    s_out=line.split()
                    n_out=len(s_out)
                    if flag==0:
                        variable_svr_mtx=np.zeros((1,n_out-8))
                        perm_force_svr=np.zeros((1,1))
                    if flag==1:
                        temp_svr=np.zeros((1,n_out-8))
                        temp_perm=np.zeros((1,1))
                        for i in range(6, n_out-1):
                            if i<n_out-2:
                                temp_svr[0,i-6]=s_out[i]
                            else:
                                temp_perm[0,0]=s_out[i]
                        variable_svr_mtx=np.row_stack((variable_svr_mtx,temp_svr))
                        perm_force_svr=np.row_stack((perm_force_svr,temp_perm))
                    if s_out[0]!='Temp' and flag==0:
                        flag=1
                        for i in range(6,n_out-1):
                            if i<n_out-2:
                                variable_svr_mtx[0,i-6]=s_out[i]
                            else:
                                perm_force_svr[0,0]=s_out[i]
            f_out.close()
            sims=len(variable_svr_mtx)           
            variable_mtx=restart_mtx
            input_files=len(restart_mtx)
            num_lines=n_re
    
    #Split Matrix
    sim_run=maxjobperrun
    if sim_run>input_files:
        sim_run=input_files
    sim_mtx=np.zeros((sim_run,num_lines))
    for i in range(sim_run):
        sim_mtx[i]=variable_mtx[0]
        variable_mtx=np.delete(variable_mtx,0,0)
    
    
    #Write Restart file
    f_r=open('Restart_next.txt','w+')
    f_r.write("%d\n" %(run_no+1))
    for i in range(len(variable_mtx)):
        for j in range(num_lines):
            if j<num_lines-1:
                f_r.write("%0.4f " %variable_mtx[i,j])
            else:
                f_r.write("%0.4f\n" %variable_mtx[i,j])
    f_r.close()
    
    #Crate file to track bigcases in run
    with open('bigcases.txt','a') as f:
        f.write('\nRun %i:\n' % run_no)

    #DSMC Simulations
    if __name__ == 'automated_multy':
        
        pool = mp.Pool(poolsize)
        
        if convergence_flag:
            for specie in species_var:
                variable_force_mtx = np.array(pool.starmap(loop_process,[( cluster, typestl, stlfile, convertionfactor, TandP, i, convergence_flag,  np.asarray(np.where(np.all(sim_mtx==i,axis=1))),sims,pathmain,MainName, specie) for i in sim_mtx]))
        else:
            for specie in species_var:
                variable_force_mtx = np.array(pool.starmap(loop_process,[( cluster, typestl, stlfile, convertionfactor, i, casesConverg, convergence_flag, np.asarray(np.where(np.all(sim_mtx==i,axis=1))),sims,pathmain,MainName, specie) for i in sim_mtx]))

        pool.close()
        
        ###### I am hereeeee!!!
        variable_force_mtx=variable_force_mtx.reshape(sim_run,-1)

        print('variable_force_mtx is', variable_force_mtx)
        flagbigcases = 1
        for i in variable_force_mtx: # 2  1  50  0.770525  100  597.692  1984  1.45948644e-09  7.90595398e-06
            if i[6] == 0:
       
                domain_extend, temp_number = i[4]*10**(-6), i[5]
                print('Entering pp_parallel_fast for temp%i'%temp_number)
                bigcasesresults = pp_parallel_fast(temp_number,domain_extend,pathmain) # return T, (P1+P2)/2, K , Perm_force, timefloat
       
                bigcasesresults = np.hstack((i[:5],bigcasesresults))
    
                path_member_log=pathmain+MainName
                member_log=os.path.join(path_member_log,'member_log.txt' )
                pathf=os.path.join(pathmain+'/Results_multi/dsmc_temp%d' %(temp_number),'log.txt')
                
                f_member=open(member_log,'a')
                with open (pathf,'w') as f_log:
                    for j in range(0,len(bigcasesresults)-1):
                        if j==0:
        
                            f_member.write(' %i    ' %temp_number)
                            f_member.write('Gas   ')
                            
                        else:
                            f_log.write('%0.3f    ' %bigcasesresults[j-1])
                            f_member.write('%0.4f    ' %bigcasesresults[j-1])
                            
                    f_log.write('%s  ' %bigcasesresults[-2])
                    f_member.write('%s    ' %bigcasesresults[-2])
                    
                    f_log.write('%s' %bigcasesresults[-1])
                    f_member.write('%s\n' %bigcasesresults[-1])
                f_member.close()
        
        # for i in range(sim_run):
        #     Perm_force[i,0]=variable_force_mtx[i,num_lines]
        # variable_force_mtx=np.delete(variable_force_mtx,-1,axis=1)
        
        # #Stack matrix for SVR
        # if run_no==1:
        #     variable_svr_mtx=np.zeros((sim_run,num_lines))
        #     perm_force_svr=np.zeros((sim_run,1))
        #     variable_svr_mtx=variable_force_mtx
        #     perm_force_svr=Perm_force
        # else:
        #     variable_svr_mtx=np.row_stack((variable_svr_mtx,variable_force_mtx))
        #     perm_force_svr=np.row_stack((perm_force_svr,Perm_force))
        
        # #Write Restart file
        # f_r=open('Restart.txt','w+')
        # f_r.write("%d\n" %(run_no+1))
        # for i in range(len(variable_mtx)):
        #     for j in range(num_lines):
        #         if j<num_lines-1:
        #             f_r.write("%0.4f " %variable_mtx[i,j])
        #         else:
        #             f_r.write("%0.4f\n" %variable_mtx[i,j])
        # f_r.close()
        
        email = 'bdde223@uky.edu'
       # phone = '8594901117@txt.att.net' # works for ATT
        alert('JobStatus', 'The simulations are done!', email)
       # alert('JobStatus', 'The simulations are done!', phone)
        
        print('All done!!!')
