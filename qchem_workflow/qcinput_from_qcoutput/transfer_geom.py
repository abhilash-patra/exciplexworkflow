from CHactive.outputs import QCOutput
from pymatgen.io.qchem.inputs import QCInput
import numpy as np
import glob

for inp in glob.glob('*.out'):
    label = inp[:-4]

    new_rems = {
    "BASIS":"6-31g*",
    "JOB_TYPE":"sp",
    "METHOD":"wB97X-D3",
    "geom_opt_max_cycles":"200",
    "MEM_STATIC":"5000",
    "MEM_TOTAL":"192000",
    "SOLVENT_METHOD":"PCM",
    "SCF_CONVERGENCE":"8",
    "SYMMETRY":"FALSE",
    "SYMMETRY_IGNORE":"1",
    "THRESH":"14"
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


    output = QCOutput(filename = inp)
    opt_geom = output.data['initial_molecule'] #output.data['molecule_from_last_geometry']
    NewInput = QCInput(molecule = opt_geom, rem = new_rems, pcm = new_pcms, solvent = new_solvents)
    NewInput.write_file(label +'.inp')
