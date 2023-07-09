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

rotation_flag1 = True #rotation not applied if False
theta1 = 83.66 #in degrees. Apply active rotation to galaxy1 by angle theta1 about y-axis in anti-clockwise direction
phi1 = -64.71 #in degrees. Apply active rotation to galaxy1 by angle phi1 about z-axis in anti-clockwise direction

#position translations

pos_translation_flag1 = True #position translation not applied if False
x1 = -4.13 #in kpc. Active translation of x-coordinates of galaxy1 by an amount x1
y1 = -4.68 #in kpc. Active translation of y-coordinates of galaxy1 by an amount y1
z1 = 0 #in kpc. Active translation of z-coordinates of galaxy1 by an amount z1

#velocity translations

vel_translation_flag1 = True #velocity translation not applied if False
v_x1 = 12.65 #in km/s. Active translation of v_x-coordinates of galaxy1 by an amount v_x1
v_y1 = 0.16 #in km/s. Active translation of v_y-coordinates of galaxy1 by an amount v_y1
v_z1 = 0.03 #in km/s. Active translation of v_z-coordinates of galaxy1 by an amount v_z1

#galaxy2

#rotation
#order of operation is first theta and then phi
#rotations are applied before translations
#rotations will be applied in the frame of the galaxy

rotation_flag2 = True #rotation not applied if False
theta2 = 90.00 #in degrees. Apply active rotation to galaxy2 by angle theta2 about y-axis in anti-clockwise direction
phi2 = 170.00 #in degrees. Apply active rotation to galaxy2 by angle phi2 about z-axis in anti-clockwise direction

#position translations
pos_translation_flag2 = True #position translation not applied if False
x2 = 35.57 #in kpc. Active translation of x-coordinates of galaxy2 by an amount x2
y2 = 40.31 #in kpc. Active translation of y-coordinates of galaxy2 by an amount y2
z2 = 0 #in kpc. Active translation of z-coordinates of galaxy2 by an amount z2

#velocity translations
vel_translation_flag2 = True #velocity translation not applied if False
v_x2 = -108.98 #in km/s. Active translation of v_x-coordinates of galaxy2 by an amount v_x2
v_y2 = -1.24 #in km/s. Active translation of v_y-coordinates of galaxy2 by an amount v_y2
v_z2 = -0.13 #in km/s. Active translation of v_z-coordinates of galaxy2 by an amount v_z2


