#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 15:06:03 2021

@author: vijaybmohan
"""

import fileinput
import sys
import time
import os
import numpy as np
import multiprocessing as mp
from DSMC_script import loop_process
from pyDOE import lhs
from alert import alert


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
if run_no==1:        
    f_input=open('textinput.txt','r')
    s=tuple()
    input_files=450 #how many points
    num_lines = sum(1 for line in open('textinput.txt'))
    input_mtx=np.zeros((num_lines,2))
    i=0
    dim_mtx=np.zeros((num_lines,1))
    for num,line in enumerate(f_input,1):
        s=line.split()
        for j in range(2):
            input_mtx[i,j]=float(s[j])
        i+=1
        if num == 1:
            bounds = [(input_mtx[num-1,0],input_mtx[num-1,1])]
        else:
            bounds.append((input_mtx[num-1,0],input_mtx[num-1,1]))
    f_input.close()

    sims=0
#Variable Matrix
    variable_mtx=np.zeros((input_files,num_lines))
    LHD = lhs(num_lines,samples=input_files,criterion="m")
    for i in range(input_files):
        for j in range(num_lines):
            variable_mtx[i,j] = input_mtx[j,0] + LHD[i,j]*(input_mtx[j,1] - input_mtx[j,0])
else:
    s_out=tuple()
    flag=0
    f_out=open('member_log.txt','r')
    for num, line in enumerate(f_out, 1):
        if num>1:
            s_out=line.split()
            n_out=len(s_out)
            if flag==0:
                variable_svr_mtx=np.zeros((1,n_out-3))
                perm_force_svr=np.zeros((1,1))
            if flag==1:
                temp_svr=np.zeros((1,n_out-3))
                temp_perm=np.zeros((1,1))
                for i in range(2,n_out):
                    if i<n_out-2:
                        temp_svr[0,i-2]=s_out[i]
                    else:
                        temp_perm[0,0]=s_out[i]
                variable_svr_mtx=np.row_stack((variable_svr_mtx,temp_svr))
                perm_force_svr=np.row_stack((perm_force_svr,temp_perm))
            if s_out[0]!='Temp' and flag==0:
                flag=1
                for i in range(2,n_out):
                    if i<n_out-2:
                        variable_svr_mtx[0,i-2]=s_out[i]
                    else:
                        perm_force_svr[0,0]=s_out[i]
    f_out.close()
    sims=len(variable_svr_mtx)           
    variable_mtx=restart_mtx
    input_files=len(restart_mtx)
    num_lines=n_re

#Split Matrix
sim_run=25 #limit per batch
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

#DSMC Simulations
variable_force_mtx=np.zeros((sim_run,num_lines+1))
Perm_force=np.zeros((sim_run,1))
if __name__ == '__main__':
    pool = mp.Pool(5)
    variable_force_mtx = np.array(pool.starmap(loop_process,[(i,np.asarray(np.where(np.all(sim_mtx==i,axis=1))),sims,pathmain,MainName) for i in sim_mtx]))
    pool.close()
variable_force_mtx=variable_force_mtx.reshape(sim_run,-1)

for i in range(sim_run):
    Perm_force[i,0]=variable_force_mtx[i,num_lines]
variable_force_mtx=np.delete(variable_force_mtx,-1,axis=1)

#Stack matrix for SVR
if run_no==1:
    variable_svr_mtx=np.zeros((sim_run,num_lines))
    perm_force_svr=np.zeros((sim_run,1))
    variable_svr_mtx=variable_force_mtx
    perm_force_svr=Perm_force
else:
    variable_svr_mtx=np.row_stack((variable_svr_mtx,variable_force_mtx))
    perm_force_svr=np.row_stack((perm_force_svr,Perm_force))


email = 'lch285@g.uky.edu'
phone = '8594901117@txt.att.net' # works for ATT
alert('JobStatus', 'The simulations are done!', email)
alert('JobStatus', 'The simulations are done!', phone)
