#Parameter file for combine_ics
#This file should be in the same directory as the executable combine_ics.py

#specifying files
ic1 = "./ics_lmc_hern_dm_halo_exp_disk.hdf5" #first ic file 
ic2 = "./ics_smc_hern_dm_halo_exp_disk.hdf5" #second ic file
ic  = "./ics_lmc_smc_hern_dm_halo_exp_disk.hdf5" #combined file

#specifying transformation parameters

#galaxy1

#rotation
#order of operation is first theta and then phi
#rotations are applied before translations
#rotations will be applied in the frame of the galaxy

theta1 = 83.66 #in degrees. Apply active rotation to galaxy1 by angle theta1 about y-axis in anti-clockwise direction
phi1 = -64.71 #in degrees. Apply active rotation to galaxy1 by angle phi1 about z-axis in anti-clockwise direction

#position translations

x1 = -4.13 #in kpc. Translate galaxy1 to new x-coordinate x1
y1 = -4.68 #in kpc. Translate galaxy1 to new y-coordinate y1
z1 = 0 #in kpc. Translate galaxy1 to new z-coordinate z1

#velocity translations
v_x1 = 12.65 #in km/s. Translate galaxy1 to new v_x-coordinate v_x1
v_y1 = 0.16 #in km/s. Translate galaxy1 to new v_y-coordinate v_y1
v_z1 = 0.03 #in km/s. Translate galaxy1 to new v_z-coordinate v_z1

#galaxy2

#rotation
#order of operation is first theta and then phi
#rotations are applied before translations
#rotations will be applied in the frame of the galaxy

theta2 = 90.00 #in degrees. Apply active rotation to galaxy2 by angle theta2 about y-axis in anti-clockwise direction
phi2 = 170.00 #in degrees. Apply active rotation to galaxy2 by angle phi2 about z-axis in anti-clockwise direction

#position translations

x2 = 35.57 #in kpc. Translate galaxy2 to new x-coordinate x2
y2 = 40.31 #in kpc. Translate galaxy2 to new y-coordinate y2
z2 = 0 #in kpc. Translate galaxy2 to new z-coordinate z2

#velocity translations
v_x2 = -108.98 #in km/s. Translate galaxy2 to new v_x-coordinate v_x2
v_y2 = -1.24 #in km/s. Translate galaxy2 to new v_y-coordinate v_y2
v_z2 = -0.13 #in km/s. Translate galaxy2 to new v_z-coordinate v_z2


