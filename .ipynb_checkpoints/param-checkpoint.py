#Parameter file for combine_ics
#This file should be in the same directory as the executable combine_ics.py

#specifying files
ic1 = "<path>" #full path of first ic file 
ic2 = "<path>" #full path of second ic file
ic  = "<path>" #full path to combined file

#particle types (mention which ones in the form of a list)
'''
Currently supported particle types:
PartType0 -> gas
PartType1 -> halo
PartType2 -> disk
PartType3 -> bulge
PartType5 -> black hole
'''
part_types1 = ['PartType0', 'PartType1', 'PartType2', 'PartType5'] #particle types for ic1
part_types2 = ['PartType0', 'PartType1', 'PartType2', 'PartType5'] #particle types for ic2

#specifying transformation parameters

#galaxy1

#rotation
#order of operation is first theta and then phi
#rotations are applied before translations
#rotations will be applied in the frame of the galaxy

rotation_flag1 = True #rotation not applied if False
theta1 = 0 #in degrees. Apply active rotation to galaxy1 by angle theta1 about y-axis in anti-clockwise direction
phi1 = 0 #in degrees. Apply active rotation to galaxy1 by angle phi1 about z-axis in anti-clockwise direction

#position translations

pos_translation_flag1 = True #position translation not applied if False
x1 = 0 #in kpc. Active translation of x-coordinates of galaxy1 by an amount x1
y1 = 0 #in kpc. Active translation of y-coordinates of galaxy1 by an amount y1
z1 = 0 #in kpc. Active translation of z-coordinates of galaxy1 by an amount z1

#velocity translations

vel_translation_flag1 = True #velocity translation not applied if False
v_x1 = 0 #in km/s. Active translation of v_x-coordinates of galaxy1 by an amount v_x1
v_y1 = 0 #in km/s. Active translation of v_y-coordinates of galaxy1 by an amount v_y1
v_z1 = 0 #in km/s. Active translation of v_z-coordinates of galaxy1 by an amount v_z1

#galaxy2

#rotation
#order of operation is first theta and then phi
#rotations are applied before translations
#rotations will be applied in the frame of the galaxy

rotation_flag2 = True #rotation not applied if False
theta2 = 0 #in degrees. Apply active rotation to galaxy2 by angle theta2 about y-axis in anti-clockwise direction
phi2 = 0 #in degrees. Apply active rotation to galaxy2 by angle phi2 about z-axis in anti-clockwise direction

#position translations
pos_translation_flag2 = True #position translation not applied if False
x2 = 0 #in kpc. Active translation of x-coordinates of galaxy2 by an amount x2
y2 = 0 #in kpc. Active translation of y-coordinates of galaxy2 by an amount y2
z2 = 0 #in kpc. Active translation of z-coordinates of galaxy2 by an amount z2

#velocity translations
vel_translation_flag2 = True #velocity translation not applied if False
v_x2 = 0 #in km/s. Active translation of v_x-coordinates of galaxy2 by an amount v_x2
v_y2 = 0 #in km/s. Active translation of v_y-coordinates of galaxy2 by an amount v_y2
v_z2 = 0 #in km/s. Active translation of v_z-coordinates of galaxy2 by an amount v_z2

#end
