from qcfw.input_eda import QCInputEDA
from pymatgen.io.qchem.outputs import QCOutput
import numpy as np



import glob
new_rems = {
"BASIS":"6-311G**",
"METHOD":"wb97xd",
"JOB_TYPE":"EDA",
"CIS_N_ROOTS":"10",
"CIS_SINGLETS":"1",
"CIS_TRIPLETS":"0",
"EX_EDA":"TRUE",
"SYMMETRY":"FALSE",
"SYM_IGNORE":"TRUE",
"THRESH":"14",
}

fragment = {
"1":"10",
}



for inp in glob.glob('e*.out'):
    label = inp[0:-4]
    output = QCOutput(filename = inp)
    opt_geom = output.data['molecule_from_last_geometry']
    NewInput = QCInputEDA(molecule = opt_geom, rem = new_rems, frgm_cis_n_roots = fragment)
    NewInput.write_file(label +'.inp')


