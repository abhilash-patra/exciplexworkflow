import os
import glob
import numpy as np
from readfile import *


rems = [
    '$rem',"\n"
'BASIS  =  GEN',"\n"
'ECP = GEN',"\n"
'METHOD = wb97xd',"\n"
'job_type = optimization',"\n"
'geom_opt_max_cycles = 700',"\n"
'max_scf_cycles = 300',"\n"
'MEM_STATIC = 5000',"\n"
'MEM_TOTAL = 192000',"\n"
'SYMMETRY  =  FALSE',"\n"
'SYMMETRY_IGNORE  =  true',"\n"
'scf_convergence = 8',"\n"
'THRESH = 14',"\n"
'UNRESTRICTED = 1',"\n"
'$end',"\n","\n"
]

new_basis = [
'$ecp',"\n"
'Cu',"\n"
'srsc',"\n"
'****',"\n"
'$end',"\n"
' ',"\n"
'$basis',"\n"
'Cu',"\n"
'srsc',"\n"
'****',"\n"
'O',"\n"
'6-311+G*',"\n"
'****',"\n"
'C',"\n"
'6-311+G*',"\n"
'****',"\n"
'N',"\n"
'6-311+G*',"\n"
'****',"\n"
'H',"\n"
'6-311+G*',"\n"
'****',"\n"
'F',"\n"
'6-311+G*',"\n"
'****',"\n"
'S',"\n"
'6-311+G*',"\n"
'****',"\n"
'$end',"\n"
]

for inp in glob.glob('*.xyz'):
    name = inp[0:-4]
    atoms, coordinates = read_xyz(inp)
    data =np.column_stack((atoms, coordinates))
    file = open(name + ".inp", "w")
    file.write("$molecule"+"\n")
    file.write("2 1"+"\n")  # Here 2 is for charge and 1 is for multiplicity. Change these values accordingly
    data1 = '\n'.join(['   '.join(['{:2}'.format(item) for item in row])for row in data])
    file.write(str(data1) + "\n")
    file.write("$end"+"\n")
    file.write("\n")
    for rem in rems:
        file.writelines(rem)

    for basis in new_basis:
        file.writelines(basis)   #This same format can be followed to add other sections in the input like solvent, fragments, cdft etc.

    file.close()
