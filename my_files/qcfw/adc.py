#!/usr/bin/env python
# coding: utf-8

# In[1]:

# Modified pymatgen code for excited state parsing for Q-Chem.5.4. General features are copied from pymatgen
#__author__ = "Abhilash Patra"


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

#-----------------------EXCITED STATE PARSING--------------------------------------------------------------------------------
#


#ADC(2)
        adc2_pcm0 = []

        temp_adc2_pcm0 = read_pattern(self.text, {"key": r"\sADC\(2\)\/ptSS-PCM\(PTE\):\s*([\d\.]+)\s*eV"}).get("key")
        self.data["adc2_pcm0"] = float(temp_adc2_pcm0[0][0])



# Open the file
#        for line in filename:
        # Try to find a match in the line
#            match = re.findall(r"ADC\(2\)\/ptSS-PCM\(PTE\):\s*([\d\.]+)\s*eV", line)
        # If a match was found
#            if match:
#               self.data["adc2_pcm0"] = float(match[0])

#-------------------Excited state EDA-------------------------------------------------------------------------------------

        ex1_eda = []
        ex2_eda = []
        ex3_eda = []
        ex4_eda = []
        ex5_eda = []


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

           dw_eda_4_temp = read_pattern(self.text, {"key": r"\sExcited state  4:\s*([\d\-\.]+)\seV\s\("}).get("key")
           dw_frz_4 = dw_eda_4_temp[0][0]
           dw_pol_4 = dw_eda_4_temp[1][0]
           dw_ct_4 = dw_eda_4_temp[2][0]
           dw_eda_4 = dw_frz_4, dw_pol_4, dw_ct_4
           ex4_eda.append(dw_eda_4)
           self.data["exc4_eda_values"] = ex4_eda

           dw_eda_5_temp = read_pattern(self.text, {"key": r"\sExcited state  5:\s*([\d\-\.]+)\seV\s\("}).get("key")
           dw_frz_5 = dw_eda_5_temp[0][0]
           dw_pol_5 = dw_eda_5_temp[1][0]
           dw_ct_5 = dw_eda_5_temp[2][0]
           dw_eda_5 = dw_frz_5, dw_pol_5, dw_ct_5
           ex4_eda.append(dw_eda_4)
           self.data["exc5_eda_values"] = ex5_eda

        except Exception: 
             pass


#        ex_eda = []
#
#
#        for i in range(1, 6):
#            try:
#               dw_eda_temp = read_pattern(self.text, {"key": r"\sExcited state  {}:\s*([\d\-\.]+)\seV\s\(".format(i)}).get("key")
#               dw_frz = dw_eda_temp[0][0]
#               dw_pol = dw_eda_temp[1][0]
#               dw_ct = dw_eda_temp[2][0]
#               dw_eda = dw_frz, dw_pol, dw_ct
#               ex_eda.append(dw_eda)
#               self.data["exc{}_eda_values".format(i)] = ex_eda
#            except Exception:
#                pass
#
 

        self._read_ex1_charges()  #Charges on atoms in excited state 1
        self._read_ex2_charges()  #Charges on atoms in excited state 2
        self._read_ex3_charges()  #Charges on atoms in excited state 3
        self._read_ex4_charges()  #Charges on atoms in excited state 4
        self._read_ex5_charges()  #Charges on atoms in excited state 5

#--------------------------------------------------------Exciton-----------------------------------------------------------

        electron_hole_distance_temp = read_pattern(self.text, {"key": r"\s\|\<r\_e \- r\_h\>\| \[Ang\]:                  *([\d\-\.]+)"}).get("key")
        self.data["electron_hole_distance"] = electron_hole_distance_temp
        
        hole_size_temp = read_pattern(self.text, {"key": r"\sHole size \[Ang\]:                      *([\d\-\.]+)"}).get("key")
        self.data["hole_size"] = hole_size_temp

        electron_size_temp = read_pattern(self.text, {"key": r"\sElectron size \[Ang\]:                  *([\d\-\.]+)"}).get("key")
        self.data["electron_size"] = electron_size_temp

        rms_separation_temp = read_pattern(self.text, {"key": r"\sRMS electron-hole separation \[Ang\]:   *([\d\-\.]+)"}).get("key")
        self.data["rms_separation"] = rms_separation_temp

        covariance_temp = read_pattern(self.text, {"key": r"\sCovariance\(r\_h\,\s* r\_e\) \[Ang\^2\]:\s*([\d\-\.]+)"}).get("key")
        self.data["covariance"] = covariance_temp

        correlation_coef_temp = read_pattern(self.text, {"key": r"\sCorrelation coefficient:              *([\d\-\.]+)"}).get("key")
        self.data["correlation_coefficient"] = correlation_coef_temp


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

#        self._read_charges()
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


#------------------------------------------Reactant_product distance----------


        temp_react_product = read_pattern(self.text, {"key": r"Total\sR\s-->\sP\sDistance\s+([\d\-\.]+)"}).get("key")
        if temp_react_product is None:
            self.data["reactant_product_dist"] = None
        else:
            self.data["reactant_product_dist"] = float(temp_react_product[0][0])


#---------------------------EX CHARGE PARSING--------------------------------------
    def _read_ex1_charges(self):
        """
        Parses excited state 1 Mulliken charges.
        """
        if self.data.get("unrestricted", []):
            header_pattern = (
                r"RPA Excited State  1:  Mulliken Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+"
                r"Spin\s\(a\.u\.\)\s+\-+"
            )
            table_pattern = r"\s+\d+\s\w+\s+([\d\-\.]+)\s+([\d\-\.]+)"
            footer_pattern = r"\s\s\-+\s+Sum of atomic charges"
        else:
            header_pattern = r"RPA Excited State  1:  Mulliken Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+\-+"
            table_pattern = r"\s+\d+\s\w+\s+([\d\-\.]+)"
            footer_pattern = r"\s\s\-+\s+Sum of atomic charges"

        temp_mulliken = read_table_pattern(self.text, header_pattern, table_pattern, footer_pattern)
        real_mulliken = []
        for one_mulliken in temp_mulliken:
            if self.data.get("unrestricted", []):
                temp = np.zeros(shape=(len(one_mulliken), 2))
                for ii, entry in enumerate(one_mulliken):
                    temp[ii, 0] = float(entry[0])
                    temp[ii, 1] = float(entry[1])
            else:
                temp = np.zeros(len(one_mulliken))
                for ii, entry in enumerate(one_mulliken):
                    temp[ii] = float(entry[0])
            real_mulliken += [temp]

        self.data["EX1_Mulliken"] = real_mulliken


    def _read_ex2_charges(self):
        """
        Parses excited state 2 Mulliken charges.
        """
        if self.data.get("unrestricted", []):
            header_pattern = (
                r"RPA Excited State  2:  Mulliken Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+"
                r"Spin\s\(a\.u\.\)\s+\-+"
            )
            table_pattern = r"\s+\d+\s\w+\s+([\d\-\.]+)\s+([\d\-\.]+)"
            footer_pattern = r"\s\s\-+\s+Sum of atomic charges"
        else:
            header_pattern = r"RPA Excited State  2:  Mulliken Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+\-+"
            table_pattern = r"\s+\d+\s\w+\s+([\d\-\.]+)"
            footer_pattern = r"\s\s\-+\s+Sum of atomic charges"


        temp_mulliken = read_table_pattern(self.text, header_pattern, table_pattern, footer_pattern)
        real_mulliken = []
        for one_mulliken in temp_mulliken:
            if self.data.get("unrestricted", []):
                temp = np.zeros(shape=(len(one_mulliken), 2))
                for ii, entry in enumerate(one_mulliken):
                    temp[ii, 0] = float(entry[0])
                    temp[ii, 1] = float(entry[1])
            else:
                temp = np.zeros(len(one_mulliken))
                for ii, entry in enumerate(one_mulliken):
                    temp[ii] = float(entry[0])
            real_mulliken += [temp]

        self.data["EX2_Mulliken"] = real_mulliken



    def _read_ex3_charges(self):
        """
        Parses excited state 3 Mulliken charges.
        """
        if self.data.get("unrestricted", []):
            header_pattern = (
                r"RPA Excited State  3:  Mulliken Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+"
                r"Spin\s\(a\.u\.\)\s+\-+"
            )
            table_pattern = r"\s+\d+\s\w+\s+([\d\-\.]+)\s+([\d\-\.]+)"
            footer_pattern = r"\s\s\-+\s+Sum of atomic charges"
        else:
            header_pattern = r"RPA Excited State  3:  Mulliken Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+\-+"
            table_pattern = r"\s+\d+\s\w+\s+([\d\-\.]+)"
            footer_pattern = r"\s\s\-+\s+Sum of atomic charges"


        temp_mulliken = read_table_pattern(self.text, header_pattern, table_pattern, footer_pattern)
        real_mulliken = []
        for one_mulliken in temp_mulliken:
            if self.data.get("unrestricted", []):
                temp = np.zeros(shape=(len(one_mulliken), 2))
                for ii, entry in enumerate(one_mulliken):
                    temp[ii, 0] = float(entry[0])
                    temp[ii, 1] = float(entry[1])
            else:
                temp = np.zeros(len(one_mulliken))
                for ii, entry in enumerate(one_mulliken):
                    temp[ii] = float(entry[0])
            real_mulliken += [temp]

        self.data["EX3_Mulliken"] = real_mulliken


    def _read_ex4_charges(self):
        """
        Parses excited state 4 Mulliken charges.
        """
        if self.data.get("unrestricted", []):
            header_pattern = (
                r"RPA Excited State  4:  Mulliken Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+"
                r"Spin\s\(a\.u\.\)\s+\-+"
            )
            table_pattern = r"\s+\d+\s\w+\s+([\d\-\.]+)\s+([\d\-\.]+)"
            footer_pattern = r"\s\s\-+\s+Sum of atomic charges"
        else:
            header_pattern = r"RPA Excited State  4:  Mulliken Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+\-+"
            table_pattern = r"\s+\d+\s\w+\s+([\d\-\.]+)"
            footer_pattern = r"\s\s\-+\s+Sum of atomic charges"


        temp_mulliken = read_table_pattern(self.text, header_pattern, table_pattern, footer_pattern)
        real_mulliken = []
        for one_mulliken in temp_mulliken:
            if self.data.get("unrestricted", []):
                temp = np.zeros(shape=(len(one_mulliken), 2))
                for ii, entry in enumerate(one_mulliken):
                    temp[ii, 0] = float(entry[0])
                    temp[ii, 1] = float(entry[1])
            else:
                temp = np.zeros(len(one_mulliken))
                for ii, entry in enumerate(one_mulliken):
                    temp[ii] = float(entry[0])
            real_mulliken += [temp]

        self.data["EX4_Mulliken"] = real_mulliken

    def _read_ex5_charges(self):
        """
        Parses excited state 5 Mulliken charges.
        """
        if self.data.get("unrestricted", []):
            header_pattern = (
                r"RPA Excited State  5:  Mulliken Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+"
                r"Spin\s\(a\.u\.\)\s+\-+"
            )
            table_pattern = r"\s+\d+\s\w+\s+([\d\-\.]+)\s+([\d\-\.]+)"
            footer_pattern = r"\s\s\-+\s+Sum of atomic charges"
        else:
            header_pattern = r"RPA Excited State  5:  Mulliken Net Atomic Charges\s+Atom\s+Charge \(a\.u\.\)\s+\-+"
            table_pattern = r"\s+\d+\s\w+\s+([\d\-\.]+)"
            footer_pattern = r"\s\s\-+\s+Sum of atomic charges"


        temp_mulliken = read_table_pattern(self.text, header_pattern, table_pattern, footer_pattern)
        real_mulliken = []
        for one_mulliken in temp_mulliken:
            if self.data.get("unrestricted", []):
                temp = np.zeros(shape=(len(one_mulliken), 2))
                for ii, entry in enumerate(one_mulliken):
                    temp[ii, 0] = float(entry[0])
                    temp[ii, 1] = float(entry[1])
            else:
                temp = np.zeros(len(one_mulliken))
                for ii, entry in enumerate(one_mulliken):
                    temp[ii] = float(entry[0])
            real_mulliken += [temp]

        self.data["EX5_Mulliken"] = real_mulliken


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

