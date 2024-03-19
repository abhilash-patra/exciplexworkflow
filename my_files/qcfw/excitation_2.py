#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re
import numpy as np
from monty.io import zopen
from monty.json import MSONable, jsanitize
from collections import defaultdict
import csv



class QCOutputES(MSONable):
    """
    Class to parse QChem output files.
    """

    def __init__(self, filename: str):
        """
        Args:
            filename (str): Filename to parse
        """
        self.filename = filename
#        print(filename)
        self.data = {}  # type: Dict[str, Any]
        self.data["errors"] = []
        self.data["warnings"] = {}
        self.text = ""
        with zopen(filename, mode="rt", encoding="ISO-8859-1") as f:
            self.text = f.read()

        # Check if output file contains multiple output files. If so, print an error message and exit
#        self.data["multiple_outputs"] = read_pattern(
#            self.text, {"key": r"Job\s+\d+\s+of\s+(\d+)\s+"}, terminate_on_match=True
#        ).get("key")
#        if self.data.get("multiple_outputs") is not None:
#            if self.data.get("multiple_outputs") != [["1"]]:
#                raise ValueError(
#                    "ERROR: multiple calculation outputs found in file "
#                    + filename
#                    + ". Please instead call QCOutput.mulitple_outputs_from_file(QCOutput,'"
#                    + filename
#                    + "')"
#                )
#
#-----------------------EXCITED STATE PARSING--------------------------------------------------------------------------------
#

        temp_cis_n_roots = read_pattern(self.text, {"key": r"CIS_N_ROOTS  =\s*([\d\-\.]+)"}).get("key")
        if temp_cis_n_roots is None:
            temp_cis_n_roots = read_pattern(self.text, {"key": r"cis_n_roots =\s*([\d\-\.]+)"}).get("key")


        self.data["cis_n_roots"] = int(temp_cis_n_roots[0][0])
        temp_cis_singlet = read_pattern(self.text, {"key": r"CIS_SINGLETS  =\s*([\d\-\.]+)"}).get("key")
        if temp_cis_singlet is None:
            temp_cis_singlet = read_pattern(self.text, {"key": r"cis_singlets = \s*([\d\-\.]+)"}).get("key")


        self.data["cis_singlet"] = int(temp_cis_singlet[0][0])
        temp_cis_triplet = read_pattern(self.text, {"key": r"CIS_TRIPLETS  =\s*([\d\-\.]+)"}).get("key")
        if temp_cis_triplet is None:
            temp_cis_triplet = read_pattern(self.text, {"key": r"cis_triplets  =\s*([\d\-\.]+)"}).get("key")



        if temp_cis_triplet is None:
            self.data["cis_triplet"] = 1
        else:
            self.data["cis_triplet"] = int(temp_cis_triplet[0][0])
        
        cis_n_roots = self.data["cis_n_roots"]
            
        excitation_energy = []
        excitation_strength = []
        excited_state_energy = []
        
        temp_strength = read_pattern(self.text, {"key": r"\sStrength   :     ([\d\-\.]+)"}).get("key")
        for jj in np.arange(1, cis_n_roots+1, 1):
            temp_excitation_energy = read_pattern(self.text, {"key": r"Excited\sstate\s+"+str(jj)+r": excitation energy \(eV\) = \s*([\d\-\.]+)"}).get("key")
            temp_excitation_energy_2 = temp_excitation_energy[-1]
            temp_excited_state_energy = read_pattern(self.text, {"key": r"Total energy for state\s+"+str(jj)+r": \s*([\d\-\.]+)\sau"}).get("key")
            temp_excited_state_energy_2 = temp_excited_state_energy[-1]
            strength = temp_strength[-(cis_n_roots+1-jj)]
            excitation_energy.insert(jj, temp_excitation_energy_2)
            excited_state_energy.insert(jj, temp_excited_state_energy_2)
            excitation_strength.insert(jj, strength)
            
        self.data["excitation_energy"] = excitation_energy
        self.data["excitation_strength"] = excitation_strength
        self.data["excited_state_energy"] = excited_state_energy 


#-------------------Excited stae EDA-------------------------------------------------------------------------------------

        ex1_eda = []
        ex2_eda = []
        ex3_eda = []



        try:
           dw_eda_1_temp = read_pattern(self.text, {"key": r"\sExcited state  1:\s*([\d\-\.]+)\seV\s\("}).get("key") 
           dw_frz_1 = dw_eda_1_temp[0][0]
           dw_pol_1 = dw_eda_1_temp[1][0]
           dw_ct_1 = dw_eda_1_temp[2][0]
           dw_eda_1 = dw_frz_1, dw_pol_1, dw_ct_1
           ex1_eda.append(dw_eda_1)
           self.data["exc1_eda_values"] = ex1_eda

           dw_eda_2_temp = read_pattern(self.text, {"key": r"\sExcited state  2:\s*([\d\-\.]+)\seV\s\("}).get("key") 
           dw_frz_2 = dw_eda_2_temp[0][0]
           dw_pol_2 = dw_eda_2_temp[1][0]
           dw_ct_2 = dw_eda_2_temp[2][0]
           dw_eda_2 = dw_frz_2, dw_pol_2, dw_ct_2
           ex2_eda.append(dw_eda_2)
           self.data["exc2_eda_values"] = ex2_eda

           dw_eda_3_temp = read_pattern(self.text, {"key": r"\sExcited state  3:\s*([\d\-\.]+)\seV\s\("}).get("key") 
           dw_frz_3 = dw_eda_3_temp[0][0]
           dw_pol_3 = dw_eda_3_temp[1][0]
           dw_ct_3 = dw_eda_3_temp[2][0]
           dw_eda_3 = dw_frz_3, dw_pol_3, dw_ct_3
           ex3_eda.append(dw_eda_3)
           self.data["exc3_eda_values"] = ex3_eda
        except Exception: 
             pass




#----------------------------------------------------------------------------------------------------------------- ---------
#---------------------------------------------------------------------------------------------------------------------------

        temp_Total_energy = read_pattern(
            self.text, {"key": r"Total energy in the final basis set =\s*([\d\-\.]+)"}
        ).get("key")
        if temp_Total_energy is not None:
            if len(temp_Total_energy) == 1:
                self.data["Total_energy_in_the_final_basis_set"] = float(temp_Total_energy[0][0])
            else:
                Total_energy = np.zeros(len(temp_Total_energy))
                for ii, val in enumerate(temp_Total_energy):
                    Total_energy[ii] = float(val[0])
                self.data["Total_energy_in_the_final_basis_set"] = Total_energy


        self._read_SCF()
        # Parse the final energy
        temp_final_energy = read_pattern(self.text, {"key": r"Final\senergy\sis\s+([\d\-\.]+)"}).get("key")
        if temp_final_energy is None:
            self.data["final_energy"] = None
        else:
            self.data["final_energy"] = float(temp_final_energy[0][0])
        
            # Check if the calculation is a single point. If so, parse the relevant output
        self.data["single_point_job"] = read_pattern(
            self.text,
            {"key": r"(?i)\s*job(?:_)*type\s*(?:=)*\s*sp"},
            terminate_on_match=True,
        ).get("key")
        if self.data.get("single_point_job", []):
            self._read_single_point_data()

            
    def _read_SCF(self):
        """
        Parses both old and new SCFs.
        """
        if self.data.get("using_GEN_SCFMAN", []):
            if "SCF_failed_to_converge" in self.data.get("errors"):
                footer_pattern = r"^\s*gen_scfman_exception: SCF failed to converge"
            else:
                footer_pattern = r"^\s*\-+\n\s+SCF time"
            header_pattern = (
                r"^\s*\-+\s+Cycle\s+Energy\s+(?:(?:DIIS)*\s+[Ee]rror)*(?:RMS Gradient)*\s+\-+"
                r"(?:\s*\-+\s+OpenMP\s+Integral\s+computing\s+Module\s+"
                r"(?:Release:\s+version\s+[\d\-\.]+\,\s+\w+\s+[\d\-\.]+\, "
                r"Q-Chem Inc\. Pittsburgh\s+)*\-+)*\n"
            )
            table_pattern = (
                r"(?:\s*Nonlocal correlation = [\d\-\.]+e[\d\-]+)*"
                r"(?:\s*Inaccurate integrated density:\n\s+Number of electrons\s+=\s+[\d\-\.]+\n\s+"
                r"Numerical integral\s+=\s+[\d\-\.]+\n\s+Relative error\s+=\s+[\d\-\.]+\s+\%\n)*\s*\d+\s+"
                r"([\d\-\.]+)\s+([\d\-\.]+)e([\d\-\.\+]+)(?:\s+Convergence criterion met)*"
                r"(?:\s+Preconditoned Steepest Descent)*(?:\s+Roothaan Step)*(?:\s+"
                r"(?:Normal\s+)*BFGS [Ss]tep)*(?:\s+LineSearch Step)*(?:\s+Line search: overstep)*"
                r"(?:\s+Dog-leg BFGS step)*(?:\s+Line search: understep)*"
                r"(?:\s+Descent step)*(?:\s+Done DIIS. Switching to GDM)*"
                r"(?:\s*\-+\s+Cycle\s+Energy\s+(?:(?:DIIS)*\s+[Ee]rror)*"
                r"(?:RMS Gradient)*\s+\-+(?:\s*\-+\s+OpenMP\s+Integral\s+computing\s+Module\s+"
                r"(?:Release:\s+version\s+[\d\-\.]+\,\s+\w+\s+[\d\-\.]+\, "
                r"Q-Chem Inc\. Pittsburgh\s+)*\-+)*\n)*"
                r"(?:(\n\s*[a-z\dA-Z_\s/]+\.C|\n\s*GDM)::WARNING energy changes are now smaller than effective "
                r"accuracy\.\s*(\n\s*[a-z\dA-Z_\s/]+\.C|\n\s*GDM)::\s+calculation will continue, but THRESH s"
                r"hould be increased\s*"
                r"(\n\s*[a-z\dA-Z_\s/]+\.C|\n\s*GDM)::\s+or SCF_CONVERGENCE decrea"
                r"sed\.\s*(\n\s*[a-z\dA-Z_\s/]+\.C|\n\s*GDM)::\s+effective_thresh = [\d\-\.]+e[\d\-]+)*"
            )
        else:
            if "SCF_failed_to_converge" in self.data.get("errors"):
                footer_pattern = r"^\s*\d+\s*[\d\-\.]+\s+[\d\-\.]+E[\d\-\.]+\s+Convergence\s+failure\n"
            else:
                footer_pattern = r"^\s*\-+\n"
            header_pattern = r"^\s*\-+\s+Cycle\s+Energy\s+DIIS Error\s+\-+\n"
            table_pattern = (
                r"(?:\s*Inaccurate integrated density:\n\s+Number of electrons\s+=\s+[\d\-\.]+\n\s+"
                r"Numerical integral\s+=\s+[\d\-\.]+\n\s+Relative error\s+=\s+[\d\-\.]+\s+\%\n)*\s*\d+\s*"
                r"([\d\-\.]+)\s+([\d\-\.]+)E([\d\-\.\+]+)(?:\s*\n\s*cpu\s+[\d\-\.]+\swall\s+[\d\-\.]+)*"
                r"(?:\nin dftxc\.C, eleTot sum is:[\d\-\.]+, tauTot is\:[\d\-\.]+)*"
                r"(?:\s+Convergence criterion met)*(?:\s+Done RCA\. Switching to DIIS)*"
                r"(?:\n\s*Warning: not using a symmetric Q)*"
                r"(?:\nRecomputing EXC\s*[\d\-\.]+\s*[\d\-\.]+\s*[\d\-\.]+"
                r"(?:\s*\nRecomputing EXC\s*[\d\-\.]+\s*[\d\-\.]+\s*[\d\-\.]+)*)*"
            )

        temp_scf = read_table_pattern(self.text, header_pattern, table_pattern, footer_pattern)
        real_scf = []
        for one_scf in temp_scf:
            temp = np.zeros(shape=(len(one_scf), 2))
            for ii, entry in enumerate(one_scf):
                temp[ii, 0] = float(entry[0])
                temp[ii, 1] = float(entry[1]) * 10 ** float(entry[2])
            real_scf += [temp]

        self.data["SCF"] = real_scf

        temp_thresh_warning = read_pattern(
            self.text,
            {
                "key": r"\n[a-zA-Z_\s/]+\.C::WARNING energy changes are now smaller than effective accuracy"
                r"\.\n[a-zA-Z_\s/]+\.C::\s+calculation will continue, but THRESH should be increased\n"
                r"[a-zA-Z_\s/]+\.C::\s+or SCF_CONVERGENCE decreased\. \n"
                r"[a-zA-Z_\s/]+\.C::\s+effective_thresh = ([\d\-\.]+e[\d\-]+)"
            },
        ).get("key")
        if temp_thresh_warning is not None:
            if len(temp_thresh_warning) == 1:
                self.data["warnings"]["thresh"] = float(temp_thresh_warning[0][0])
            else:
                thresh_warning = np.zeros(len(temp_thresh_warning))
                for ii, entry in enumerate(temp_thresh_warning):
                    thresh_warning[ii] = float(entry[0])
                self.data["warnings"]["thresh"] = thresh_warning

        temp_SCF_energy = read_pattern(self.text, {"key": r"SCF   energy in the final basis set =\s*([\d\-\.]+)"}).get(
            "key"
        )
        if temp_SCF_energy is not None:
            if len(temp_SCF_energy) == 1:
                self.data["SCF_energy_in_the_final_basis_set"] = float(temp_SCF_energy[0][0])
            else:
                SCF_energy = np.zeros(len(temp_SCF_energy))
                for ii, val in enumerate(temp_SCF_energy):
                    SCF_energy[ii] = float(val[0])
                self.data["SCF_energy_in_the_final_basis_set"] = SCF_energy

        temp_Total_energy = read_pattern(
            self.text, {"key": r"Total energy in the final basis set =\s*([\d\-\.]+)"}
        ).get("key")
        if temp_Total_energy is not None:
            if len(temp_Total_energy) == 1:
                self.data["Total_energy_in_the_final_basis_set"] = float(temp_Total_energy[0][0])
            else:
                Total_energy = np.zeros(len(temp_Total_energy))
                for ii, val in enumerate(temp_Total_energy):
                    Total_energy[ii] = float(val[0])
                self.data["Total_energy_in_the_final_basis_set"] = Total_energy     
                
                
    def _read_single_point_data(self):
        """
        Parses final free energy information from single-point calculations.
        """
        temp_dict = read_pattern(
            self.text,
            {"final_energy": r"\s*Total\s+energy in the final basis set\s+=\s*([\d\-\.]+)"},
        )

        if temp_dict.get("final_energy") is None:
            self.data["final_energy"] = None
        else:
            # -1 in case of pcm
            # Two lines will match the above; we want final calculation
            self.data["final_energy"] = float(temp_dict.get("final_energy")[-1][0])
            
            
def read_pattern(text_str, patterns, terminate_on_match=False, postprocess=str):

    compiled = {key: re.compile(pattern, re.MULTILINE | re.DOTALL) for key, pattern in patterns.items()}
    matches = defaultdict(list)
    for key, pattern in compiled.items():
        for match in pattern.finditer(text_str):
            matches[key].append([postprocess(i) for i in match.groups()])
            if terminate_on_match:
                break
    return matches


def read_table_pattern(
    text_str,
    header_pattern,
    row_pattern,
    footer_pattern,
    postprocess=str,
    attribute_name=None,
    last_one_only=False,
):
    table_pattern_text = header_pattern + r"\s*(?P<table_body>(?:" + row_pattern + r")+)\s*" + footer_pattern
    table_pattern = re.compile(table_pattern_text, re.MULTILINE | re.DOTALL)
    rp = re.compile(row_pattern)
    data = {}
    tables = []
    for mt in table_pattern.finditer(text_str):
        table_body_text = mt.group("table_body")
        table_contents = []
        for ml in rp.finditer(table_body_text):
            d = ml.groupdict()
            if len(d) > 0:
                processed_line = {k: postprocess(v) for k, v in d.items()}
            else:
                processed_line = [postprocess(v) for v in ml.groups()]
            table_contents.append(processed_line)
        tables.append(table_contents)
    if last_one_only:
        retained_data = tables[-1]
    else:
        retained_data = tables
    if attribute_name is not None:
        data[attribute_name] = retained_data
        return data
    return retained_data

