from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import QCOutput

import glob


new_rems = {
"BASIS":"6-311g**",
"METHOD":"wb97xd",
"JOB_TYPE":"SP",
"CIS_N_ROOTS":"10",
"CIS_SINGLETS":"1",
"CIS_TRIPLETS":"0",
"RPA":"2",
"MAX_SCF_CYCLES":"200",
"SOLVENT_METHOD":"PCM",
"SYMMETRY":"FALSE",
"SYMMETRY_IGNORE":"1",
"NTO_PAIRS":"2",
"GUI":"2",
"CIS_AMPL_ANAL":"1",
"STATE_ANALYSIS":"1"
}

new_pcms = {
"HEAVYPOINTS":"590",
"METHOD":"SWIG",
"RADII":"BONDI",
"SOLVER":"INVERSION",
"THEORY":"CPCM"
}

new_solvents = {
"DIELECTRIC":"2.03",
"OPTICALDIELECTRIC":"2.01"
}



for inp in glob.glob('*.out'):
    label = inp[0:-4]
    output = QCOutput(filename = inp)
    opt_geom = output.data['molecule_from_last_geometry']
    NewInput = QCInput(molecule = opt_geom, rem = new_rems, pcm = new_pcms, solvent = new_solvents)
    NewInput.write_file(label +'.inp')

