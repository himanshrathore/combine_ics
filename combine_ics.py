#Code to combine two HDF5 initial conditions
#Author: Himansh Rathore, February 2023
#Last updated: Apr 21, 2023 by Himansh Rathore

import numpy as np
import h5py
import os
import param #parameter file param.py should be in the same directory as the script combine_ics.py

#############################function definitions#######################################################################

def translate_gal(state, frame): #active translation
    #state should be a 3 * n_particles matrix containing coordinates (position/velocity etc.) of the star particles
    #frame should be a (3, 1) array containing the coordinates (position/velocity) of the new reference frame

    return np.subtract(state.T, frame.reshape(1,3)).T

def Ry(beta): #active rotation about y-axis
    return np.array([[np.cos(beta), 0, np.sin(beta)],
                    [0, 1, 0],
                    [-np.sin(beta), 0, np.cos(beta)]])

def Rz(gamma): #active rotation about z-axis
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

N_part_types = 6 #total number of gadget particle types
print('Assuming 6 particle types')

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
ic1_part1_count = int(head1.attrs['NumPart_Total'][1])
ic1_part2_count = int(head1.attrs['NumPart_Total'][2])
ic1_part5_count = int(head1.attrs['NumPart_Total'][5])

ic2_part1_count = int(head2.attrs['NumPart_Total'][1])
ic2_part2_count = int(head2.attrs['NumPart_Total'][2])
ic2_part5_count = int(head2.attrs['NumPart_Total'][5])

N_ic1 = int(np.sum(head1.attrs['NumPart_Total'])) #total no. of particles in ic1
N_ic2 = int(np.sum(head2.attrs['NumPart_Total'])) #total no. of particles in ic2

#Assuming no gas

#creating PartType1 (DM halo)

#creating ParticleIDs
ic1_part1_pids = np.array(f1['PartType1']['ParticleIDs'], dtype = int)
ic2_part1_pids = np.array(f2['PartType1']['ParticleIDs'], dtype = int)

part1_pids = np.append(ic1_part1_pids, ic2_part1_pids + N_ic1).astype(int) 

dset = f.create_dataset('/PartType1/ParticleIDs', shape = part1_pids.shape, dtype = part1_pids.dtype, data = part1_pids)

print("Creating particle IDs of PartType1...")

#creating coordinates

ic1_part1_coord = np.array(f1['PartType1']['Coordinates'], dtype = np.float64) #halo coordinates of ic1
#performing rotation
ic1_part1_coord = np.matmul(Rz(phi1), np.matmul(Ry(theta1), ic1_part1_coord.T)).T
#performing translation
ic1_part1_coord = translate_gal(ic1_part1_coord.T, np.array([x1, y1, z1])).T

ic2_part1_coord = np.array(f2['PartType1']['Coordinates'], dtype = np.float64) #halo coordinates of ic2
#performing rotation
ic2_part1_coord = np.matmul(Rz(phi2), np.matmul(Ry(theta2), ic2_part1_coord.T)).T
#performing translation
ic2_part1_coord = translate_gal(ic2_part1_coord.T, np.array([x2, y2, z2])).T

part1_coord = np.vstack((ic1_part1_coord, ic2_part1_coord)) #combined halo coordinates

dset = f.create_dataset('/PartType1/Coordinates', shape = part1_coord.shape, dtype = part1_coord.dtype, data = part1_coord)

print("Creating coordinates of PartType1...")

#creating velocities

ic1_part1_vel = np.array(f1['PartType1']['Velocities'], dtype = np.float64) #halo velocities of ic1
#performing rotation
ic1_part1_vel = np.matmul(Rz(phi1), np.matmul(Ry(theta1), ic1_part1_vel.T)).T
#performing translation
ic1_part1_vel = translate_gal(ic1_part1_vel.T, np.array([v_x1, v_y1, v_z1])).T

ic2_part1_vel = np.array(f2['PartType1']['Velocities'], dtype = np.float64) #halo velocities of ic2
#performing rotation
ic2_part1_vel = np.matmul(Rz(phi2), np.matmul(Ry(theta2), ic2_part1_vel.T)).T
#performing translation
ic2_part1_vel = translate_gal(ic2_part1_vel.T, np.array([v_x2, v_y2, v_z2])).T

part1_vel = np.vstack((ic1_part1_vel, ic2_part1_vel)) #combined halo velocities

dset = f.create_dataset('/PartType1/Velocities', shape = part1_vel.shape, dtype = part1_vel.dtype, data = part1_vel)

print("Creating velocities of PartType1...")

#creating masses

#try getting mass from the mass list of the particle type
try:
    ic1_part1_mass = np.array(f1['PartType1']['Masses'], dtype = np.float64, ndmin=1)
    #if this does not exist, particles of this type have identical mass
    #the mass is listed in the header of the table
except KeyError:
    ic1_part1_mass = np.full(ic1_part1_count, head1.attrs['MassTable'][1]).astype(float)
    
try:
    ic2_part1_mass = np.array(f2['PartType1']['Masses'], dtype = np.float64, ndmin=1)
    #if this does not exist, particles of this type have identical mass
    #the mass is listed in the header of the table
except KeyError:
    ic2_part1_mass = np.full(ic2_part1_count, head2.attrs['MassTable'][1]).astype(float)
    
part1_masses = np.append(ic1_part1_mass, ic2_part1_mass).astype(float)

dset = f.create_dataset('/PartType1/Masses', shape = part1_masses.shape, dtype = part1_masses.dtype, data = part1_masses)

print("Creating masses of PartType1...")
    
#creating PartType2 (stellar disk)

#creating ParticleIDs

#creating ParticleIDs
ic1_part2_pids = np.array(f1['PartType2']['ParticleIDs'], dtype = int)
ic2_part2_pids = np.array(f2['PartType2']['ParticleIDs'], dtype = int)

part2_pids = np.append(ic1_part2_pids, ic2_part2_pids + N_ic1).astype(int) 

dset = f.create_dataset('/PartType2/ParticleIDs', shape = part2_pids.shape, dtype = part2_pids.dtype, data = part2_pids)

print("Creating particle IDs of PartType2...")

#creating coordinates

ic1_part2_coord = np.array(f1['PartType2']['Coordinates'], dtype = np.float64) #disk coordinates of ic1
#performing rotation
ic1_part2_coord = np.matmul(Rz(phi1), np.matmul(Ry(theta1), ic1_part2_coord.T)).T
#performing translation
ic1_part2_coord = translate_gal(ic1_part2_coord.T, np.array([x1, y1, z1])).T

ic2_part2_coord = np.array(f2['PartType2']['Coordinates'], dtype = np.float64) #disk coordinates of ic2
#performing rotation
ic2_part2_coord = np.matmul(Rz(phi2), np.matmul(Ry(theta2), ic2_part2_coord.T)).T
#performing translation
ic2_part2_coord = translate_gal(ic2_part2_coord.T, np.array([x2, y2, z2])).T

part2_coord = np.vstack((ic1_part2_coord, ic2_part2_coord)) #combined disk coordinates

dset = f.create_dataset('/PartType2/Coordinates', shape = part2_coord.shape, dtype = part2_coord.dtype, data = part2_coord)

print("Creating coordinates of PartType2...")

#creating velocities

ic1_part2_vel = np.array(f1['PartType2']['Velocities'], dtype = np.float64) #disk velocities of ic1
#performing rotation
ic1_part2_vel = np.matmul(Rz(phi1), np.matmul(Ry(theta1), ic1_part2_vel.T)).T
#performing translation
ic1_part2_vel = translate_gal(ic1_part2_vel.T, np.array([v_x1, v_y1, v_z1])).T

ic2_part2_vel = np.array(f2['PartType2']['Velocities'], dtype = np.float64) #disk velocities of ic2
#performing rotation
ic2_part2_vel = np.matmul(Rz(phi2), np.matmul(Ry(theta2), ic2_part2_vel.T)).T
#performing translation
ic2_part2_vel = translate_gal(ic2_part2_vel.T, np.array([v_x2, v_y2, v_z2])).T

part2_vel = np.vstack((ic1_part2_vel, ic2_part2_vel)) #combined halo velocities

dset = f.create_dataset('/PartType2/Velocities', shape = part2_vel.shape, dtype = part2_vel.dtype, data = part2_vel)

print("Creating velocities of PartType2...")

#creating masses

#try getting mass from the mass list of the particle type
try:
    ic1_part2_mass = np.array(f1['PartType2']['Masses'], dtype = np.float64, ndmin=1)
    #if this does not exist, particles of this type have identical mass
    #the mass is listed in the header of the table
except KeyError:
    ic1_part2_mass = np.full(ic1_part2_count, head1.attrs['MassTable'][2]).astype(float)
    
try:
    ic2_part2_mass = np.array(f2['PartType2']['Masses'], dtype = np.float64, ndmin=1)
    #if this does not exist, particles of this type have identical mass
    #the mass is listed in the header of the table
except KeyError:
    ic2_part2_mass = np.full(ic2_part2_count, head2.attrs['MassTable'][2]).astype(float)
    
part2_masses = np.append(ic1_part2_mass, ic2_part2_mass).astype(float)

dset = f.create_dataset('/PartType2/Masses', shape = part2_masses.shape, dtype = part2_masses.dtype, data = part2_masses)

print("Creating masses of PartType2...")

#creating PartType5 (central BHs)

#creating ParticleIDs
ic1_part5_pids = np.array(f1['PartType5']['ParticleIDs'], dtype = int)
ic2_part5_pids = np.array(f2['PartType5']['ParticleIDs'], dtype = int)

part5_pids = np.append(ic1_part5_pids, ic2_part5_pids + N_ic1).astype(int) 

dset = f.create_dataset('/PartType5/ParticleIDs', shape = part5_pids.shape, dtype = part5_pids.dtype, data = part5_pids)

print("Creating particle IDs of PartType5...")

#creating coordinates

ic1_part5_coord = np.array(f1['PartType5']['Coordinates'], dtype = np.float64) #BH coordinates of ic1
#performing rotation
ic1_part5_coord = np.matmul(Rz(phi1), np.matmul(Ry(theta1), ic1_part5_coord.T)).T
#performing translation
ic1_part5_coord = translate_gal(ic1_part5_coord.T, np.array([x1, y1, z1])).T

ic2_part5_coord = np.array(f2['PartType5']['Coordinates'], dtype = np.float64) #BH coordinates of ic2
#performing rotation
ic2_part5_coord = np.matmul(Rz(phi2), np.matmul(Ry(theta2), ic2_part5_coord.T)).T
#performing translation
ic2_part5_coord = translate_gal(ic2_part5_coord.T, np.array([x2, y2, z2])).T

part5_coord = np.vstack((ic1_part5_coord, ic2_part5_coord)) #combined disk coordinates

dset = f.create_dataset('/PartType5/Coordinates', shape = part5_coord.shape, dtype = part5_coord.dtype, data = part5_coord)

print("Creating coordinates of PartType5...")

#creating velocities

ic1_part5_vel = np.array(f1['PartType5']['Velocities'], dtype = np.float64) #BH velocities of ic1
#performing rotation
ic1_part5_vel = np.matmul(Rz(phi1), np.matmul(Ry(theta1), ic1_part5_vel.T)).T
#performing translation
ic1_part5_vel = translate_gal(ic1_part5_vel.T, np.array([v_x1, v_y1, v_z1])).T

ic2_part5_vel = np.array(f2['PartType5']['Velocities'], dtype = np.float64) #BH velocities of ic2
#performing rotation
ic2_part5_vel = np.matmul(Rz(phi2), np.matmul(Ry(theta2), ic2_part5_vel.T)).T
#performing translation
ic2_part5_vel = translate_gal(ic2_part5_vel.T, np.array([v_x2, v_y2, v_z2])).T

part5_vel = np.vstack((ic1_part5_vel, ic2_part5_vel)) #combined halo velocities

dset = f.create_dataset('/PartType5/Velocities', shape = part5_vel.shape, dtype = part5_vel.dtype, data = part5_vel)

print("Creating velocities of PartType5...")

#creating masses

#try getting mass from the mass list of the particle type
try:
    ic1_part5_mass = np.array(f1['PartType5']['Masses'], dtype = np.float64, ndmin=1)
    #if this does not exist, particles of this type have identical mass
    #the mass is listed in the header of the table
except KeyError:
    ic1_part5_mass = np.full(ic1_part5_count, head1.attrs['MassTable'][5]).astype(float)
    
try:
    ic2_part5_mass = np.array(f2['PartType5']['Masses'], dtype = np.float64, ndmin=1)
    #if this does not exist, particles of this type have identical mass
    #the mass is listed in the header of the table
except KeyError:
    ic2_part5_mass = np.full(ic2_part5_count, head2.attrs['MassTable'][5]).astype(float)
    
part5_masses = np.append(ic1_part5_mass, ic2_part5_mass).astype(float)

dset = f.create_dataset('/PartType5/Masses', shape = part5_masses.shape, dtype = part5_masses.dtype, data = part5_masses)

print("Creating masses of PartType5...")

print("Done.")














