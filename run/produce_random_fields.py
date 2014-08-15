#!/usr/bin/python
import numpy as np
import os
import subprocess
import glob
import shutil

#USER DEFINED PARAMETERS
N_m = 20
obs_nside = 16
Nreals = 20
lmax_array = np.array([20,50,100,200,500])
n_cores = 4
#-----------------------
n_lmax = len(lmax_array)

#Generate directories to store random field/other housekeeping commands
folder_contents = os.listdir('.')
for i in range(n_lmax):
    #Assume directories do not exist
    folder_existence = False
    #Check if directories already exist
    for j in folder_contents:
        if j=='lmax_%i'%lmax_array[i]:
            folder_existence = True
    #If no folder exists, create one
    if folder_existence == False:
        os.mkdir('lmax_%i' %lmax_array[i])
    #If folder already exists, continue loop
    else:
        continue

#Seed increment depends on Nside used
if obs_nside == 1:
    Npix = 12*obs_nside**2
    seed_incr = (710-Npix)*4*N_m

if obs_nside == 2:
    Npix = 12*obs_nside**2
    seed_incr = (2934-Npix)*4*N_m

if obs_nside == 4:
    Npix = 12*obs_nside**2
    seed_incr = (11790-Npix)*4*N_m

if obs_nside == 8:
    Npix = 12*obs_nside**2
    seed_incr = (47316-Npix)*4*N_m

if obs_nside == 16:
    Npix = 12*obs_nside**2
    seed_incr = (189276-Npix)*4*N_m

if obs_nside == 32:
    Npix = 12*obs_nside**2
    seed_incr = (757472-Npix)*4*N_m

if obs_nside == 64:
    Npix = 12*obs_nside**2
    seed_incr = (3030270-Npix)*4*N_m

if obs_nside == 128:
    Npix = 12*obs_nside**2
    seed_incr = (12118494-Npix)*4*N_m

if obs_nside == 128:
    Npix = 12*obs_nside**2
    seed_incr = (48469568-Npix)*4*N_m

seed = N_m*4 #Start seed index depends on number of modes used
seed_array = np.zeros(Nreals,dtype=np.int)
temp_seed = seed
for i in range(Nreals):
    seed_array[i]= temp_seed
    temp_seed = temp_seed + seed_incr

#First element in this array should be starting seed after origin seeds
seed_array[0] = seed

#Generate parameter file list
pfile =[]
for i in range(19):
    pfile.append(0)

#Assign parameter entries to list
pfile[0] = 'B_field_type = 7\n'
pfile[1] = 'z_diameter = 8\n'
pfile[2] = 'B_field_do_random = F\n'
pfile[3] = 'max_radius = 25\n'
pfile[4] = 'vec_size_R = 150\n'
pfile[5] = 'obs_shell_index_numb = 1\n'
pfile[6] = 'total_shell_numb = 1\n'
pfile[7] = 'do_sync_emission = false\n'
pfile[8] = 'list_use_lonlat_range=T\n'
pfile[9] = 'lon_min=-1\n'
pfile[10] = 'lon_max=360\n'
pfile[11] = 'lat_min=-90\n'
pfile[12] = 'lat_max=90\n'
pfile[13] = 'list_out_file_name = destroy_me\n'
pfile[14] = '', #start_seed entry
pfile[15] = '', #lmax entry
pfile[16] = 'obs_NSIDE = %i\n' %obs_nside #obs_NSIDE entry
pfile[17] = 'N_m = %i\n' %N_m #N_m entry
pfile[18] = '' #file_name_num entry

core_ticker=1
for i in range(n_lmax):
    for j in range(Nreals):
        pfile[14] = 'start_seed = %i\n' %seed_array[j]
        pfile[15] = 'lmax = %.2f\n' %lmax_array[i]
        pfile[18] = 'file_name_num = %i' %(j+1)
        pfile_out = open('rand_pfile%i.txt'%core_ticker,'w')
        pfile_out.writelines(pfile)
        pfile_out.close()
        core_ticker = core_ticker+1

#Begin producing random field realizations
#Check if number of realizations is divisible by number of cores
remainder = Nreals % n_cores
#Write a shell script that runs the jobs, depending on # of cores
#(Need to do this because I don't understand multiprocessing.)
if remainder==0:
    sh_script = []
    for i in range(16):
        sh_script.append(0)
    sh_script[0] = '#!/bin/bash\n'
    sh_script[1] = '\n'
    sh_script[2] = 'for a in %s\n' %np.str(lmax_array)[1:-1]
    sh_script[3] = 'do\n'
    sh_script[4] = 'let ticker=1\n'
    sh_script[5] = 'for b in `seq 1 %i`;\n' %(Nreals/n_cores)
    sh_script[6] = 'do\n'
    sh_script[7] = 'for c in `seq 1 %i`;\n' %n_cores
    sh_script[8] = 'do\n'
    sh_script[9] = './hammurabi rand_pfile${ticker}.txt &\n'
    sh_script[10] = 'let ticker=ticker+1\n'
    sh_script[11] = 'done\n'
    sh_script[12] = 'wait\n'
    sh_script[13] = 'done\n'
    sh_script[14] = 'mv field_out_* lmax_${a}\n'
    sh_script[15] = 'done\n'
if remainder>0:
    sh_script = []
    for i in range(22):
        sh_script.append(0)
    sh_script[0] = '#!/bin/bash\n'
    sh_script[1] = '\n'
    sh_script[2] = 'for a in %s\n' %np.str(lmax_array)[1:-1]
    sh_script[3] = 'do\n'
    sh_script[4] = 'let ticker=1\n'
    sh_script[5] = 'for b in `seq 1 %i`;\n' %(int(Nreals/n_cores))
    sh_script[6] = 'do\n'
    sh_script[7] = 'for c in `seq 1 %i`;\n' %n_cores
    sh_script[8] = 'do\n'
    sh_script[9] = './hammurabi rand_pfile${ticker}.txt &\n'
    sh_script[10] = 'let ticker=ticker+1\n'
    sh_script[11] = 'done\n'
    sh_script[12] = 'wait\n'
    sh_script[13] = 'done\n'
    sh_script[14] = 'for d in `seq 1 %i`;\n' %(Nreals % n_cores)
    sh_script[15] = 'do\n'
    sh_script[16] = './hammurabi rand_pfile${ticker}.txt &\n'
    sh_script[17] = 'let ticker=ticker+1\n'
    sh_script[18] = 'done\n'
    sh_script[19] = 'wait\n'
    sh_script[20] = 'mv field_out_* lmax_${a}\n'
    sh_script[21] = 'done\n'

shscr_out = open('run_job.sh','w')
shscr_out.writelines(sh_script)
shscr_out.close()
subprocess.call(['chmod','700','run_job.sh'])
subprocess.call('./run_job.sh')

#---------- KRF NORMALIZATION SCRIPT ----------#
#::::: Sean Quinn
#::::: July 23, 2014
#------------------ USAGE ------------------#
#Requires random field output from Hammurabi
#Format of this should be: Bx  By  Bz  x  y  z
#where the coordinates have units of kpc
#When used at the terminal this reads in the file
#which should have a format like: field_out_X.txt
#where X is an integer. The script then removes duplicate
#origin entries and computes the normalization factor 
#to give the realization Brms=1
#The properly normalized field is saved to a file: nfield_out_X.txt
#------------------------------------------#

norm_consts = np.zeros((Nreals,2))
for i in range(n_lmax):
    os.chdir('lmax_%i' %lmax_array[i])
    for j in range(Nreals):
        bx,by,bz,x,y,z = np.loadtxt("field_out_%i.txt" %(j+1),unpack=True)
        #Determine the non-origin indices of the loaded vectors
        diff_array = (x+8.5)+(y+0)+(z+0) #These elements are zero iff. x,y,z is the sun coordinate
        non_orig_indices = np.where(np.abs(diff_array) > 0.)
        #Create vectors of all points excluding origin
        bx2 = np.empty(0)
        by2 = np.empty(0)
        bz2 = np.empty(0)
        bx2 = np.ravel(np.take(bx,non_orig_indices))
        by2 = np.ravel(np.take(by,non_orig_indices))
        bz2 = np.ravel(np.take(bz,non_orig_indices))
        #Now add 1 copy of origin values. This entire procedure ensures
        #origin value is not multiply counted
        bx2 = np.append(bx2,bx[0])
        by2 = np.append(by2,by[0])
        bz2 = np.append(bz2,bz[0])
        #Now compute the rms of the field magnitude, and store it
        brms = np.sqrt(np.sum(bx2**2+by2**2+bz2**2)/len(bx2))
        norm_consts[j,0] = j+1
        norm_consts[j,1] = brms
        #Now normalize the field grid to have Brms=1
        Bnorm = np.zeros((len(bx),3))
        Bnorm[:,0] = bx / brms
        Bnorm[:,1] = by / brms
        Bnorm[:,2] = bz / brms
        #Save the normalized grid to a separate file. The file name is
        #the same with n prepended.
        np.savetxt("nfield_out_%i.txt" %(j+1), Bnorm, delimiter='\t', fmt="%.15f")
    os.chdir('..')



