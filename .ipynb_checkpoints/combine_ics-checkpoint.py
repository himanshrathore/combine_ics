#Code to combine two Gadget-4 based HDF5 initial conditions
#Author: Himansh Rathore, February 2023
#Last updated: Aug 14, 2025 by Himansh Rathore

import numpy as np
import h5py
import os
import param #parameter file param.py should be in the same directory as the script combine_ics.py

#############################function definitions#######################################################################

def translate_gal(state, frame): #active translation
    #state should be a 3 * n_particles matrix containing coordinates (position/velocity etc.) of the star particles
    #frame should be a (3, 1) array containing the coordinates (position/velocity) of the new reference frame

    return np.subtract(state.T, frame.reshape(1,3)).T

def Ry(beta): #active rotation matrix about y-axis
    return np.array([[np.cos(beta), 0, np.sin(beta)],
                    [0, 1, 0],
                    [-np.sin(beta), 0, np.cos(beta)]])

def Rz(gamma): #active rotation matrix about z-axis
    return np.array([[np.cos(gamma), -np.sin(gamma), 0],
                    [np.sin(gamma), np.cos(gamma), 0],
                    [0, 0, 1]])

########################################################################################################################

print('Starting code...')

#Accessing the files
ic1 = param.ic1 #ic1
ic2 = param.ic2 #ic2
ic = param.ic #combined ic

#Reading the IC files
f1 = h5py.File(ic1, 'r') #reading ic1
print('Reading galaxy1: ', ic1)
f2 = h5py.File(ic2, 'r') #reading ic2
print('Reading galaxy2: ', ic2)

#Creating the combined IC
#check if the file already exists
if(os.path.exists(ic)):
    #removing file
    os.remove(ic)
    print('Combined ic file already exists. Will remove it and create a new one.')

f = h5py.File(ic, 'w')
print('Starting to create combined IC: ', ic)

#Reading particle types
part_types1 = param.part_types1
part_types2 = param.part_types2
print('Reading particle types...')

#Reading transformation parameters
theta1 = param.theta1*np.pi/180
phi1 = param.phi1*np.pi/180
x1 = -param.x1
y1 = -param.y1
z1 = -param.z1
v_x1 = -param.v_x1
v_y1 = -param.v_y1
v_z1 = -param.v_z1
theta2 = param.theta2*np.pi/180
phi2 = param.phi2*np.pi/180
x2 = -param.x2
y2 = -param.y2
z2 = -param.z2
v_x2 = -param.v_x2
v_y2 = -param.v_y2
v_z2 = -param.v_z2
print('Reading transformation parameters...')

#Reading flags
rotation_flag1 = param.rotation_flag1
rotation_flag2 = param.rotation_flag2
pos_translation_flag1 = param.pos_translation_flag1
pos_translation_flag2 = param.pos_translation_flag2
vel_translation_flag1 = param.vel_translation_flag1
vel_translation_flag2 = param.vel_translation_flag2
print('Reading transformation flags...')

N_part_types = 6 #total number of gadget particle types
print('Assuming a total of 6 possible particle types')

#header of ic1
head1 = f1['Header']
#header of ic2
head2 = f2['Header']

#header of combined file
head = f.create_group('/Header')

head.attrs['BoxSize'] = head1.attrs['BoxSize'] #keeping same as ic1
head.attrs['MassTable'] = np.zeros(N_part_types) #since particles of a given type can have different masses
head.attrs['NumFilesPerSnapshot'] = head1.attrs['NumFilesPerSnapshot'] #keeping same as ic1
head.attrs['NumPart_ThisFile'] = head1.attrs['NumPart_ThisFile'] + head2.attrs['NumPart_ThisFile']
head.attrs['NumPart_Total'] = head1.attrs['NumPart_Total'] + head2.attrs['NumPart_Total']
head.attrs['Redshift'] = head1.attrs['Redshift'] #keeping same as ic1
head.attrs['Time'] = 0.0 #since this is the ic

print('Creating header of combined ic')

#keeping number count of different particle types

ic1_part0_count = int(head1.attrs['NumPart_Total'][0])
ic1_part1_count = int(head1.attrs['NumPart_Total'][1])
ic1_part2_count = int(head1.attrs['NumPart_Total'][2])
ic1_part3_count = int(head1.attrs['NumPart_Total'][3])
ic1_part5_count = int(head1.attrs['NumPart_Total'][5])

ic2_part0_count = int(head2.attrs['NumPart_Total'][0])
ic2_part1_count = int(head2.attrs['NumPart_Total'][1])
ic2_part2_count = int(head2.attrs['NumPart_Total'][2])
ic2_part3_count = int(head2.attrs['NumPart_Total'][3])
ic2_part5_count = int(head2.attrs['NumPart_Total'][5])

N_ic1 = int(np.sum(head1.attrs['NumPart_Total'])) #total no. of particles in ic1
N_ic2 = int(np.sum(head2.attrs['NumPart_Total'])) #total no. of particles in ic2

####################################################################################################################################################

#creating PartType0 (gas)

if(('PartType0' in part_types1) or ('PartType0' in part_types2)): #do nothing if gas is not present in either IC's

    #creating ParticleIDs

    if('PartType0' in part_types1):
        ic1_part0_pids = np.array(f1['PartType0']['ParticleIDs'], dtype = int)
    else:
        ic1_part0_pids = np.array([])

    if('PartType0' in part_types2):
        ic2_part0_pids = np.array(f2['PartType0']['ParticleIDs'], dtype = int)
    else:
        ic2_part0_pids = np.array([])

    part0_pids = np.append(ic1_part0_pids, ic2_part0_pids + N_ic1).astype(int) 

    dset = f.create_dataset('/PartType0/ParticleIDs', shape = part0_pids.shape, dtype = part0_pids.dtype, data = part0_pids)

    print("Creating particle IDs of PartType0...")

    #creating coordinates
    
    if('PartType0' in part_types1):
    
        ic1_part0_coord = np.array(f1['PartType0']['Coordinates'], dtype = np.float64) #gas coordinates of ic1
        #performing rotation
        if(rotation_flag1 == True):
            ic1_part0_coord = np.matmul(Rz(phi1), np.matmul(Ry(theta1), ic1_part0_coord.T)).T
        #performing translation
        if(pos_translation_flag1 == True):
            ic1_part0_coord = translate_gal(ic1_part0_coord.T, np.array([x1, y1, z1])).T
    
    else:
        
        ic1_part0_coord = np.array([])
        
    if('PartType0' in part_types2):

        ic2_part0_coord = np.array(f2['PartType0']['Coordinates'], dtype = np.float64) #gas coordinates of ic2
        #performing rotation
        if(rotation_flag2 == True):
            ic2_part0_coord = np.matmul(Rz(phi2), np.matmul(Ry(theta2), ic2_part0_coord.T)).T
        #performing translation
        if(pos_translation_flag2 == True):
            ic2_part0_coord = translate_gal(ic2_part0_coord.T, np.array([x2, y2, z2])).T
            
    else:
        
        ic2_part0_coord = np.array([])

    part0_coord = np.vstack((ic1_part0_coord, ic2_part0_coord)) #combined gas coordinates

    dset = f.create_dataset('/PartType0/Coordinates', shape = part0_coord.shape, dtype = part0_coord.dtype, data = part0_coord)

    print("Creating coordinates of PartType0...")

    #creating velocities
    
    if('PartType0' in part_types1):

        ic1_part0_vel = np.array(f1['PartType0']['Velocities'], dtype = np.float64) #gas velocities of ic1
        #performing rotation
        if(rotation_flag1 == True):
            ic1_part0_vel = np.matmul(Rz(phi1), np.matmul(Ry(theta1), ic1_part0_vel.T)).T
        #performing translation
        if(vel_translation_flag1 == True):
            ic1_part0_vel = translate_gal(ic1_part0_vel.T, np.array([v_x1, v_y1, v_z1])).T
            
    else:
        
        ic1_part0_vel = np.array([])
        
    if('PartType0' in part_types2):

        ic2_part0_vel = np.array(f2['PartType0']['Velocities'], dtype = np.float64) #gas velocities of ic2
        #performing rotation
        if(rotation_flag2 == True):
            ic2_part0_vel = np.matmul(Rz(phi2), np.matmul(Ry(theta2), ic2_part0_vel.T)).T
        #performing translation
        if(vel_translation_flag2 == True):
            ic2_part0_vel = translate_gal(ic2_part0_vel.T, np.array([v_x2, v_y2, v_z2])).T
            
    else:
        
        ic2_part0_vel = np.array([])

    part0_vel = np.vstack((ic1_part0_vel, ic2_part0_vel)) #combined gas velocities

    dset = f.create_dataset('/PartType0/Velocities', shape = part0_vel.shape, dtype = part0_vel.dtype, data = part0_vel)

    print("Creating velocities of PartType0...")

    #creating masses
    
    if('PartType0' in part_types1):
    
        #try getting mass from the mass list of the particle type
        try:
            ic1_part0_mass = np.array(f1['PartType0']['Masses'], dtype = np.float64, ndmin=1)
            #if this does not exist, particles of this type have identical mass
            #the mass is listed in the header of the table
        except KeyError:
            ic1_part0_mass = np.full(ic1_part0_count, head1.attrs['MassTable'][0]).astype(float)
            
    else:
        
        ic1_part0_mass = np.array([])
        
    if('PartType0' in part_types2):

        try:
            ic2_part0_mass = np.array(f2['PartType0']['Masses'], dtype = np.float64, ndmin=1)
            #if this does not exist, particles of this type have identical mass
            #the mass is listed in the header of the table
        except KeyError:
            ic2_part0_mass = np.full(ic2_part0_count, head2.attrs['MassTable'][0]).astype(float)
            
    else:
        
        ic2_part0_mass = np.array([])

    part0_masses = np.append(ic1_part0_mass, ic2_part0_mass).astype(float)

    dset = f.create_dataset('/PartType0/Masses', shape = part0_masses.shape, dtype = part0_masses.dtype, data = part0_masses)

    print("Creating masses of PartType0...")

####################################################################################################################################################

#creating PartType1 (DM halo)

if(('PartType1' in part_types1) or ('PartType1' in part_types2)): #do nothing if dm is not present in either IC's

    #creating ParticleIDs

    if('PartType1' in part_types1):
        ic1_part1_pids = np.array(f1['PartType1']['ParticleIDs'], dtype = int)
    else:
        ic1_part1_pids = np.array([])

    if('PartType1' in part_types2):
        ic2_part1_pids = np.array(f2['PartType1']['ParticleIDs'], dtype = int)
    else:
        ic2_part1_pids = np.array([])

    part1_pids = np.append(ic1_part1_pids, ic2_part1_pids + N_ic1).astype(int) 

    dset = f.create_dataset('/PartType1/ParticleIDs', shape = part1_pids.shape, dtype = part1_pids.dtype, data = part1_pids)

    print("Creating particle IDs of PartType1...")

    #creating coordinates

    if('PartType1' in part_types1):

        ic1_part1_coord = np.array(f1['PartType1']['Coordinates'], dtype = np.float64) #halo coordinates of ic1
        #performing rotation
        if(rotation_flag1 == True):
            ic1_part1_coord = np.matmul(Rz(phi1), np.matmul(Ry(theta1), ic1_part1_coord.T)).T
        #performing translation
        if(pos_translation_flag1 == True):
            ic1_part1_coord = translate_gal(ic1_part1_coord.T, np.array([x1, y1, z1])).T

    else:

        ic1_part1_coord = np.array([])

    if('PartType1' in part_types2):

        ic2_part1_coord = np.array(f2['PartType1']['Coordinates'], dtype = np.float64) #halo coordinates of ic2
        #performing rotation
        if(rotation_flag2 == True):
            ic2_part1_coord = np.matmul(Rz(phi2), np.matmul(Ry(theta2), ic2_part1_coord.T)).T
        #performing translation
        if(pos_translation_flag2 == True):
            ic2_part1_coord = translate_gal(ic2_part1_coord.T, np.array([x2, y2, z2])).T

    else:

        ic2_part1_coord = np.array([])

    part1_coord = np.vstack((ic1_part1_coord, ic2_part1_coord)) #combined halo coordinates

    dset = f.create_dataset('/PartType1/Coordinates', shape = part1_coord.shape, dtype = part1_coord.dtype, data = part1_coord)

    print("Creating coordinates of PartType1...")

    #creating velocities
    
    if('PartType1' in part_types1):

        ic1_part1_vel = np.array(f1['PartType1']['Velocities'], dtype = np.float64) #halo velocities of ic1
        #performing rotation
        if(rotation_flag1 == True):
            ic1_part1_vel = np.matmul(Rz(phi1), np.matmul(Ry(theta1), ic1_part1_vel.T)).T
        #performing translation
        if(vel_translation_flag1 == True):
            ic1_part1_vel = translate_gal(ic1_part1_vel.T, np.array([v_x1, v_y1, v_z1])).T
            
    else:
        
        ic1_part1_vel = np.array([])
        
    if('PartType1' in part_types2):

        ic2_part1_vel = np.array(f2['PartType1']['Velocities'], dtype = np.float64) #halo velocities of ic2
        #performing rotation
        if(rotation_flag2 == True):
            ic2_part1_vel = np.matmul(Rz(phi2), np.matmul(Ry(theta2), ic2_part1_vel.T)).T
        #performing translation
        if(vel_translation_flag2 == True):
            ic2_part1_vel = translate_gal(ic2_part1_vel.T, np.array([v_x2, v_y2, v_z2])).T
            
    else:
        
        ic2_part1_vel = np.array([])

    part1_vel = np.vstack((ic1_part1_vel, ic2_part1_vel)) #combined halo velocities

    dset = f.create_dataset('/PartType1/Velocities', shape = part1_vel.shape, dtype = part1_vel.dtype, data = part1_vel)

    print("Creating velocities of PartType1...")

    #creating masses
    
    if('PartType1' in part_types1):

        #try getting mass from the mass list of the particle type
        try:
            ic1_part1_mass = np.array(f1['PartType1']['Masses'], dtype = np.float64, ndmin=1)
            #if this does not exist, particles of this type have identical mass
            #the mass is listed in the header of the table
        except KeyError:
            ic1_part1_mass = np.full(ic1_part1_count, head1.attrs['MassTable'][1]).astype(float)
            
    else:
        
        ic1_part1_mass = np.array([])
        
    if('PartType1' in part_types2):

        try:
            ic2_part1_mass = np.array(f2['PartType1']['Masses'], dtype = np.float64, ndmin=1)
            #if this does not exist, particles of this type have identical mass
            #the mass is listed in the header of the table
        except KeyError:
            ic2_part1_mass = np.full(ic2_part1_count, head2.attrs['MassTable'][1]).astype(float)
            
    else:
        
        ic2_part1_mass = np.array([])

    part1_masses = np.append(ic1_part1_mass, ic2_part1_mass).astype(float)

    dset = f.create_dataset('/PartType1/Masses', shape = part1_masses.shape, dtype = part1_masses.dtype, data = part1_masses)

    print("Creating masses of PartType1...")

####################################################################################################################################################

#creating PartType2 (stellar disk)

if(('PartType2' in part_types1) or ('PartType2' in part_types2)): #do nothing if disk is not present in either IC's

    #creating ParticleIDs
    if('PartType2' in part_types1):
        ic1_part2_pids = np.array(f1['PartType2']['ParticleIDs'], dtype = int)
    else:
        ic1_part2_pids = np.array([])
    if('PartType2' in part_types2):
        ic2_part2_pids = np.array(f2['PartType2']['ParticleIDs'], dtype = int)
    else:
        ic2_part2_pids = np.array([])

    part2_pids = np.append(ic1_part2_pids, ic2_part2_pids + N_ic1).astype(int) 

    dset = f.create_dataset('/PartType2/ParticleIDs', shape = part2_pids.shape, dtype = part2_pids.dtype, data = part2_pids)

    print("Creating particle IDs of PartType2...")

    #creating coordinates

    if('PartType2' in part_types1):
        ic1_part2_coord = np.array(f1['PartType2']['Coordinates'], dtype = np.float64) #disk coordinates of ic1
        #performing rotation
        if(rotation_flag1 == True):
            ic1_part2_coord = np.matmul(Rz(phi1), np.matmul(Ry(theta1), ic1_part2_coord.T)).T
        #performing translation
        if(pos_translation_flag1 == True):
            ic1_part2_coord = translate_gal(ic1_part2_coord.T, np.array([x1, y1, z1])).T
    else:
        ic1_part2_coord = np.array([])
        
    if('PartType2' in part_types2):
        ic2_part2_coord = np.array(f2['PartType2']['Coordinates'], dtype = np.float64) #disk coordinates of ic2
        #performing rotation
        if(rotation_flag2 == True):
            ic2_part2_coord = np.matmul(Rz(phi2), np.matmul(Ry(theta2), ic2_part2_coord.T)).T
        #performing translation
        if(pos_translation_flag2 == True):
            ic2_part2_coord = translate_gal(ic2_part2_coord.T, np.array([x2, y2, z2])).T
    else:
        ic2_part2_coord = np.array([])

    part2_coord = np.vstack((ic1_part2_coord, ic2_part2_coord)) #combined disk coordinates

    dset = f.create_dataset('/PartType2/Coordinates', shape = part2_coord.shape, dtype = part2_coord.dtype, data = part2_coord)

    print("Creating coordinates of PartType2...")

    #creating velocities
    
    if('PartType2' in part_types1):
        ic1_part2_vel = np.array(f1['PartType2']['Velocities'], dtype = np.float64) #disk velocities of ic1
        #performing rotation
        if(rotation_flag1 == True):
            ic1_part2_vel = np.matmul(Rz(phi1), np.matmul(Ry(theta1), ic1_part2_vel.T)).T
        #performing translation
        if(vel_translation_flag1 == True):
            ic1_part2_vel = translate_gal(ic1_part2_vel.T, np.array([v_x1, v_y1, v_z1])).T
    else:
        ic1_part2_vel = np.array([])
        
    if('PartType2' in part_types2):
        ic2_part2_vel = np.array(f2['PartType2']['Velocities'], dtype = np.float64) #disk velocities of ic2
        #performing rotation
        if(rotation_flag2 == True):
            ic2_part2_vel = np.matmul(Rz(phi2), np.matmul(Ry(theta2), ic2_part2_vel.T)).T
        #performing translation
        if(vel_translation_flag2 == True):
            ic2_part2_vel = translate_gal(ic2_part2_vel.T, np.array([v_x2, v_y2, v_z2])).T
    else:
        ic2_part2_vel = np.array([])

    part2_vel = np.vstack((ic1_part2_vel, ic2_part2_vel)) #combined disk velocities

    dset = f.create_dataset('/PartType2/Velocities', shape = part2_vel.shape, dtype = part2_vel.dtype, data = part2_vel)

    print("Creating velocities of PartType2...")

    #creating masses

    #try getting mass from the mass list of the particle type
    if('PartType2' in part_types1):
        try:
            ic1_part2_mass = np.array(f1['PartType2']['Masses'], dtype = np.float64, ndmin=1)
            #if this does not exist, particles of this type have identical mass
            #the mass is listed in the header of the table
        except KeyError:
            ic1_part2_mass = np.full(ic1_part2_count, head1.attrs['MassTable'][2]).astype(float)
    else:
        ic1_part2_mass = np.array([])
        
    if('PartType2' in part_types2):
        try:
            ic2_part2_mass = np.array(f2['PartType2']['Masses'], dtype = np.float64, ndmin=1)
            #if this does not exist, particles of this type have identical mass
            #the mass is listed in the header of the table
        except KeyError:
            ic2_part2_mass = np.full(ic2_part2_count, head2.attrs['MassTable'][2]).astype(float)
    else:
        ic2_part2_mass = np.array([])

    part2_masses = np.append(ic1_part2_mass, ic2_part2_mass).astype(float)

    dset = f.create_dataset('/PartType2/Masses', shape = part2_masses.shape, dtype = part2_masses.dtype, data = part2_masses)

    print("Creating masses of PartType2...")

####################################################################################################################################################

#creating PartType3 (bulge)

if(('PartType3' in part_types1) or ('PartType3' in part_types2)): #do nothing if bulge is not present in either IC's

    #creating ParticleIDs
    if('PartType3' in part_types1):
        ic1_part3_pids = np.array(f1['PartType3']['ParticleIDs'], dtype = int)
    else:
        ic1_part3_pids = np.array([])
    if('PartType3' in part_types2):
        ic2_part3_pids = np.array(f2['PartType3']['ParticleIDs'], dtype = int)
    else:
        ic2_part3_pids = np.array([])

    part3_pids = np.append(ic1_part3_pids, ic2_part3_pids + N_ic1).astype(int) 

    dset = f.create_dataset('/PartType3/ParticleIDs', shape = part3_pids.shape, dtype = part3_pids.dtype, data = part3_pids)

    print("Creating particle IDs of PartType3...")

    #creating coordinates

    if('PartType3' in part_types1):
        ic1_part3_coord = np.array(f1['PartType3']['Coordinates'], dtype = np.float64) #bulge coordinates of ic1
        #performing rotation
        if(rotation_flag1 == True):
            ic1_part3_coord = np.matmul(Rz(phi1), np.matmul(Ry(theta1), ic1_part3_coord.T)).T
        #performing translation
        if(pos_translation_flag1 == True):
            ic1_part3_coord = translate_gal(ic1_part3_coord.T, np.array([x1, y1, z1])).T
    else:
        ic1_part3_coord = np.array([])
        
    if('PartType3' in part_types2):
        ic2_part3_coord = np.array(f2['PartType3']['Coordinates'], dtype = np.float64) #bulge coordinates of ic2
        #performing rotation
        if(rotation_flag2 == True):
            ic2_part3_coord = np.matmul(Rz(phi2), np.matmul(Ry(theta2), ic2_part3_coord.T)).T
        #performing translation
        if(pos_translation_flag2 == True):
            ic2_part3_coord = translate_gal(ic2_part3_coord.T, np.array([x2, y2, z2])).T
    else:
        ic2_part3_coord = np.array([])

    part3_coord = np.vstack((ic1_part3_coord, ic2_part3_coord)) #combined bulge coordinates

    dset = f.create_dataset('/PartType3/Coordinates', shape = part3_coord.shape, dtype = part3_coord.dtype, data = part3_coord)

    print("Creating coordinates of PartType3...")

    #creating velocities
    
    if('PartType3' in part_types1):
        ic1_part3_vel = np.array(f1['PartType3']['Velocities'], dtype = np.float64) #bulge velocities of ic1
        #performing rotation
        if(rotation_flag1 == True):
            ic1_part3_vel = np.matmul(Rz(phi1), np.matmul(Ry(theta1), ic1_part3_vel.T)).T
        #performing translation
        if(vel_translation_flag1 == True):
            ic1_part3_vel = translate_gal(ic1_part3_vel.T, np.array([v_x1, v_y1, v_z1])).T
    else:
        ic1_part3_vel = np.array([])
        
    if('PartType3' in part_types2):
        ic2_part3_vel = np.array(f2['PartType3']['Velocities'], dtype = np.float64) #bulge velocities of ic2
        #performing rotation
        if(rotation_flag2 == True):
            ic2_part3_vel = np.matmul(Rz(phi2), np.matmul(Ry(theta2), ic2_part3_vel.T)).T
        #performing translation
        if(vel_translation_flag2 == True):
            ic2_part3_vel = translate_gal(ic2_part3_vel.T, np.array([v_x2, v_y2, v_z2])).T
    else:
        ic2_part3_vel = np.array([])

    part3_vel = np.vstack((ic1_part3_vel, ic2_part3_vel)) #combined bulge velocities

    dset = f.create_dataset('/PartType3/Velocities', shape = part3_vel.shape, dtype = part3_vel.dtype, data = part3_vel)

    print("Creating velocities of PartType3...")

    #creating masses

    #try getting mass from the mass list of the particle type
    if('PartType3' in part_types1):
        try:
            ic1_part3_mass = np.array(f1['PartType3']['Masses'], dtype = np.float64, ndmin=1)
            #if this does not exist, particles of this type have identical mass
            #the mass is listed in the header of the table
        except KeyError:
            ic1_part3_mass = np.full(ic1_part3_count, head1.attrs['MassTable'][3]).astype(float)
    else:
        ic1_part3_mass = np.array([])
        
    if('PartType3' in part_types2):
        try:
            ic2_part3_mass = np.array(f2['PartType3']['Masses'], dtype = np.float64, ndmin=1)
            #if this does not exist, particles of this type have identical mass
            #the mass is listed in the header of the table
        except KeyError:
            ic2_part3_mass = np.full(ic2_part3_count, head2.attrs['MassTable'][3]).astype(float)
    else:
        ic2_part3_mass = np.array([])

    part3_masses = np.append(ic1_part3_mass, ic2_part3_mass).astype(float)

    dset = f.create_dataset('/PartType3/Masses', shape = part3_masses.shape, dtype = part3_masses.dtype, data = part3_masses)

    print("Creating masses of PartType3...")

####################################################################################################################################################

#creating PartType5 (central BHs)

if(('PartType5' in part_types1) or ('PartType5' in part_types2)): #do nothing if bh is not present in either IC's

    #creating ParticleIDs
    if('PartType5' in part_types1):
        ic1_part5_pids = np.array(f1['PartType5']['ParticleIDs'], dtype = int)
    else:
        ic1_part5_pids = np.array([])
    if('PartType5' in part_types2):
        ic2_part5_pids = np.array(f2['PartType5']['ParticleIDs'], dtype = int)
    else:
        ic2_part5_pids = np.array([])

    part5_pids = np.append(ic1_part5_pids, ic2_part5_pids + N_ic1).astype(int) 

    dset = f.create_dataset('/PartType5/ParticleIDs', shape = part5_pids.shape, dtype = part5_pids.dtype, data = part5_pids)

    print("Creating particle IDs of PartType5...")

    #creating coordinates
    
    if('PartType5' in part_types1):
        ic1_part5_coord = np.array(f1['PartType5']['Coordinates'], dtype = np.float64) #BH coordinates of ic1
        #performing rotation
        if(rotation_flag1 == True):
            ic1_part5_coord = np.matmul(Rz(phi1), np.matmul(Ry(theta1), ic1_part5_coord.T)).T
        #performing translation
        if(pos_translation_flag1 == True):
            ic1_part5_coord = translate_gal(ic1_part5_coord.T, np.array([x1, y1, z1])).T
    else:
        ic1_part5_coord = np.array([])
        
    if('PartType5' in part_types2):
        ic2_part5_coord = np.array(f2['PartType5']['Coordinates'], dtype = np.float64) #BH coordinates of ic2
        #performing rotation
        if(rotation_flag2 == True):
            ic2_part5_coord = np.matmul(Rz(phi2), np.matmul(Ry(theta2), ic2_part5_coord.T)).T
        #performing translation
        if(pos_translation_flag2 == True):
            ic2_part5_coord = translate_gal(ic2_part5_coord.T, np.array([x2, y2, z2])).T
    else:
        ic2_part5_coord = np.array([])

    part5_coord = np.vstack((ic1_part5_coord, ic2_part5_coord)) #combined disk coordinates

    dset = f.create_dataset('/PartType5/Coordinates', shape = part5_coord.shape, dtype = part5_coord.dtype, data = part5_coord)

    print("Creating coordinates of PartType5...")

    #creating velocities
    
    if('PartType5' in part_types1):
        ic1_part5_vel = np.array(f1['PartType5']['Velocities'], dtype = np.float64) #BH velocities of ic1
        #performing rotation
        if(rotation_flag1 == True):
            ic1_part5_vel = np.matmul(Rz(phi1), np.matmul(Ry(theta1), ic1_part5_vel.T)).T
        #performing translation
        if(vel_translation_flag1 == True):
            ic1_part5_vel = translate_gal(ic1_part5_vel.T, np.array([v_x1, v_y1, v_z1])).T
    else:
        ic1_part5_vel = np.array([])
        
    if('PartType5' in part_types2):
        ic2_part5_vel = np.array(f2['PartType5']['Velocities'], dtype = np.float64) #BH velocities of ic2
        #performing rotation
        if(rotation_flag2 == True):
            ic2_part5_vel = np.matmul(Rz(phi2), np.matmul(Ry(theta2), ic2_part5_vel.T)).T
        #performing translation
        if(vel_translation_flag2 == True):
            ic2_part5_vel = translate_gal(ic2_part5_vel.T, np.array([v_x2, v_y2, v_z2])).T
    else:
        ic2_part5_vel = np.array([])

    part5_vel = np.vstack((ic1_part5_vel, ic2_part5_vel)) #combined halo velocities

    dset = f.create_dataset('/PartType5/Velocities', shape = part5_vel.shape, dtype = part5_vel.dtype, data = part5_vel)

    print("Creating velocities of PartType5...")

    #creating masses
    
    if('PartType5' in part_types1):
    #try getting mass from the mass list of the particle type
        try:
            ic1_part5_mass = np.array(f1['PartType5']['Masses'], dtype = np.float64, ndmin=1)
            #if this does not exist, particles of this type have identical mass
            #the mass is listed in the header of the table
        except KeyError:
            ic1_part5_mass = np.full(ic1_part5_count, head1.attrs['MassTable'][5]).astype(float)
    else:
        ic1_part5_mass = np.array([])

    if('PartType5' in part_types2):
        try:
            ic2_part5_mass = np.array(f2['PartType5']['Masses'], dtype = np.float64, ndmin=1)
            #if this does not exist, particles of this type have identical mass
            #the mass is listed in the header of the table
        except KeyError:
            ic2_part5_mass = np.full(ic2_part5_count, head2.attrs['MassTable'][5]).astype(float)
    else:
        ic2_part5_mass = np.array([])
        
    part5_masses = np.append(ic1_part5_mass, ic2_part5_mass).astype(float)

    dset = f.create_dataset('/PartType5/Masses', shape = part5_masses.shape, dtype = part5_masses.dtype, data = part5_masses)

    print("Creating masses of PartType5...")

print("Done.")

#end
