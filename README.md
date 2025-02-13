# combine_ics
Combine N-body Gadget4 based HDF5 galaxy initial condition files for galaxy interaction studies.

# Author: 
Himansh Rathore, February 2023

# Usage:
The parameter file should be in the same directory as the executable "combine_ics.py" and should have the name "param.py". Then, directly execute "combine_ics.py" using python3.

Example:
<br>
python3 ./combine_ics

# Parameter file:
See the supplied "param.py" for instructions.

#Supported particle types:
Currently, the following gadget particles types are supported:
PartType0 -> gas
PartType1 -> dark matter
PartType2 -> stellar disk
PartType5 -> black hole  
