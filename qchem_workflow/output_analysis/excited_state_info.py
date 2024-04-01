from qcfw.excitation import QCOutputES
import pandas as pd
import numpy as np
import glob
from matplotlib import pyplot as plt
import csv
import os


final_energy_data = []
excitation_energy_data = []
excited_state_energy_data = []
excited_state_energy_emi_data = []
excitation_strength_1_data = []
excitation_strength_2_data = []
name_data = []

ct_data = []
frz_data = []
pol_data = []

charge_opp = []
charge_tea = []

emission_energy_data = []
oscillator_strength = []
emission_wave_length_in_nm_data = []


electron_hole_distance_data = []
hole_size_data = []
electron_size_data = []
rms_separation_data = []
covariance_data = []
correlation_coefficient_data = []

for inp in glob.glob('*.out'):
    label = inp[0:-4]
    name_data.append(label)

    os.chdir("/project/ssharada_52/patraa/qchem_workflow/output_analysis/unique_geometry")   #path for optimized geometry
    output = QCOutputES(filename = str(inp))
    final_energy = round(output.data['final_energy'],4)
    final_energy_data.append(final_energy)
    excitation_energy = str(output.data['excitation_energy'][0])
    excitation_energy = round(float(excitation_energy.replace("['", '').replace("']", '')),4)
    excitation_energy_data.append(excitation_energy)
    excited_state_energy = str(output.data['excited_state_energy'][0])
    excited_state_energy = round(float(excited_state_energy.replace("['", '').replace("']", '')),4)
    excited_state_energy_data.append(excited_state_energy)
    excitation_strength_1 = str(output.data['excitation_strength'][0])
    excitation_strength_1 = round(float(excitation_strength_1.replace("['", '').replace("']", '')),3)
    excitation_strength_1_data.append(excitation_strength_1)
    excitation_strength_2 = str(output.data['excitation_strength'][1])
    excitation_strength_2 = round(float(excitation_strength_2.replace("['", '').replace("']", '')),3)
    excitation_strength_2_data.append(excitation_strength_2)


    os.chdir("/project/ssharada_52/patraa/qchem_workflow/output_analysis/output_eda")        #path for EDA output
    output1 = QCOutputES(filename = str(inp))
    eda_values = output1.data['exc1_eda_values']
    ct_data.append(eda_values[0][2])
    frz_data.append(eda_values[0][0])
    pol_data.append(eda_values[0][1])
    

    os.chdir("/project/ssharada_52/patraa/qchem_workflow/output_analysis/charge_out")        #path for Mulliken output
    output2 = QCOutputES(filename = str(inp))
    sum0 = 0
    sum1 = 0
    for i in np.arange(0,22):
        sum0 = sum0 + output2.data['EX1_Mulliken'][-1][i]

    charge_tea.append(str(round(sum0,3)))
    for i in np.arange(22,54):
        sum1 = sum1 + output2.data['EX1_Mulliken'][-1][i]

    charge_opp.append(str(round(sum1,3)))


    os.chdir("/project/ssharada_52/patraa/qchem_workflow/output_analysis/emission_output_cyclohexane")      #path for emission calculation
    output3 = QCOutputES(filename = str(inp))   
    final_energy = round(output3.data['Total_energy_in_the_final_basis_set'][1],4)
    excited_state_energy_emi = str(output3.data['excited_state_energy'][0])
    excited_state_energy_emi = round(float(excited_state_energy_emi.replace("['", '').replace("']", '')),4)
    excited_state_energy_emi_data.append(excited_state_energy_emi)
    oscillator_strength_1 = str(output3.data['excitation_strength'][0])
    oscillator_strength_1 = oscillator_strength_1.replace("['", '').replace("']", '')
    emission_energy = round(27.2113246*(excited_state_energy_emi-final_energy),3)
    emission_wave_length_in_nm = round(45.56337117/(excited_state_energy_emi-final_energy))
    oscillator_strength.append(round(float(oscillator_strength_1),3))
    emission_energy_data.append(emission_energy)
    emission_wave_length_in_nm_data.append(emission_wave_length_in_nm)



    os.chdir("/project/ssharada_52/patraa/qchem_workflow/output_analysis/exciton_out")                       #path for exciton output
    output4 = QCOutputES(filename = str(inp))
    data_keys = ['electron_hole_distance','hole_size','electron_size','rms_separation','covariance','correlation_coefficient']
    for key in data_keys:
        value = str(output4.data[key][0])
        value = value.replace("['", '').replace("']", '')
        value = round(float(value), 4)
        globals()[f"{key}_data"].append(value)


os.chdir("/project/ssharada_52/patraa/qchem_workflow/output_analysis")
df = {'name':name_data, 'final_energy':final_energy_data, 'excited_state_energy':excited_state_energy_data, 'excitation_energy':excitation_energy_data,'excitation_strength_1':excitation_strength_1_data,'excitation_strength_2':excitation_strength_2_data,'dw_frz(eV)':frz_data, 'dw_pol(eV)':pol_data, 'dw_ct(eV)':ct_data, 'Charge_on_OPP3*':charge_opp, 'Charge_on_TEA':charge_tea,'excited_state_energy_emission':excited_state_energy_emi_data, 'emission_energy_eV':emission_energy_data, 'emission_wave_length_in_nm':emission_wave_length_in_nm_data,'Oscillator_strength':oscillator_strength, 'e_h_distance':electron_hole_distance_data}
excitation_info = pd.DataFrame(df)
excitation_info1 = excitation_info.sort_values(by='emission_energy_eV')
excitation_info1.to_csv('excited_state_data_cyclohexane.csv')
print(excitation_info1)

