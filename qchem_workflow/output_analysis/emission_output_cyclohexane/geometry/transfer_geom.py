from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import QCOutput

import glob

#convert output into QCOutput object

new_rems = {
"BASIS":"6-311G**",
"METHOD":"wb97xd",
"JOB_TYPE":"SP",
"CIS_N_ROOTS":"10",
"CIS_SINGLETS":"1",
"CIS_TRIPLETS":"0",
"RPA":"2",
"SOLVENT_METHOD":"PCM",
"SYMMETRY":"FALSE",
"SYMMETRY_IGNORE":"1",
"CIS_DYNAMIC_MEM":"TRUE",
"CIS_RELAXED_DENSITY":"TRUE",
"USE_NEW_FUNCTIONAL":"TRUE",
"PCM_PRINT":"1"
}

new_pcms = {
"ChargeSeparation":"Excited",
"StateSpecific":"1"
}

new_solvents = {
"DIELECTRIC":"2.03",
"OPTICALDIELECTRIC":"2.01"
}



solvents = [
    '$solvent',"\n"
'DIELECTRIC  2.03',"\n"
'OPTICALDIELECTRIC  2.01',"\n"
'$end',"\n"
]



rem2 = [
    '$rem',"\n"
'BASIS  =  6-311G**',"\n"
'METHOD = wb97xd',"\n"
'SCF_GUESS = READ',"\n"
'SCF_CONVERGENCE = 8',"\n"
'SOLVENT_METHOD  =  PCM',"\n"
'SYMMETRY  =  FALSE',"\n"
'SYMMETRY_IGNORE  =  1',"\n"
'SOLVENT_METHOD =  PCM',"\n"
'PCM_PRINT =  1',"\n"
'$end',"\n","\n"
]


pcm2 = [
    '$pcm',"\n"
'StateSpecific      Marcus',"\n"
'$end',"\n","\n"
]



for inp in glob.glob('e*.out'):
    label = inp[0:-4]
    output = QCOutput(filename = inp)
    opt_geom = output.data['molecule_from_last_geometry']
    NewInput = QCInput(molecule = opt_geom, rem = new_rems, pcm = new_pcms, solvent = new_solvents)
    NewInput.write_file(label +'.inp')
    file = open(label +'.inp', "a")
    file.write("\n") 
    file.write("@@@"+"\n")
    file.write("\n") 
    file.write("$molecule"+"\n")
    file.write("   READ"+"\n")
    file.write("$end"+"\n")
    file.write("\n")
    for rem in rem2:
        file.writelines(rem)
                                 
                                 
    for pcm in pcm2:
        file.writelines(pcm)
                                 
    for solvent in solvents:
        file.writelines(solvent)
    
    file.close()


