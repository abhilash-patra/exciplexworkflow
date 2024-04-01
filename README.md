# exciplexworkflow
Excited state TDDFT workflow automation using Fireworks, Custodian, and Pymatgen.

We will use MongoDB for our database. As it is not recommended to do any computation in the login node, we will use the interactive way to use one compute node for our setup and database connection. After logging in to the HPC, use the following command to access one compute node for the interactive session.

*salloc --time=02:00:00 --cpus-per-task=4  —-partition=#name of a partition*

Change the values and partition according to your needs. It will assign one compute node.

### Setting up MongoDB and other Python packages (Using conda/mamba environment)

- *module purge*                      

- *module load conda*

- *mamba init bash*

- *source ~/.bashrc*

- *mamba create --name fwworkflow*  

- *mamba activate fwworkflow*       

fwworkflow is the name of the environment. The last line activates the created environment. To install any Python package in the conda environment, search for the package on the anaconda webpage https://anaconda.org/

The following commands install the required packages.

- *mamba install conda-forge::mongodb*

- *mamba install conda-forge::fireworks python=3.10*

- *mamba install conda-forge::pymatgen*

- *mamba install conda-forge::custodian*

### Pymatgen Parsers for Q-Chem input file generation and excited state output extraction 

I have kept all the developed parsers for Q-Chem and necessary files for FireWorks in the folder ‘my_files’. Copy this file to your home and add your pythonpath to this directory as follows,

- *export PYTHONPATH=“/path-to/my_files:$PYTHONPATH"*

Make sure to replace /path-to/my_files with the actual path to your my_files directory. If you're unsure about the path, the pwd command will display it when you're inside the my_files directory. Remember to save any changes you make to .bashrc and then source the file or restart the terminal session for the changes to take effect.
- *source ~/.bashrc*

### Q-Chem Workflow Automation

### Input file creation

1. From .xyz files:
   - *cd qchem_workflow/qcinput_from_xyz*
   - *python xyz_to_inp..py*

2. From previous output files:
   - *qchem_workflow/qcinput_from_qcoutput/*
   - *transfer_geom.py*

### Output parsing and Analysis

Directory:- qchem_workflow/output_analysis

In each directory inside, I place a ‘geometry’ folder containing one example output file and a Python script named ‘transfer_geom.py’. This script is designed to create the input file.

1. Mulliken population
   - *output_analysis/charge_out/geometry*
3. Emission
   - *output_analysis/emission_output_cyclohexane/geometry*
5. Excited state EDA
   - *output_analysis/output_eda/geometry*
7. Exciton
   - *output_analysis/exciton_out/geometry*
9. Convert output to .xyz
    - *output_analysis/xyz_files*
   
Each directory (output_analysis/*) contains two sample output files. The Python script ‘excited_state_info.py’, located in the parent directory, is designed to extract all necessary information from the output files in each subdirectory and create a consolidated data file. To execute this process correctly, the paths for each directory must be updated in the Python script.





### Pymatgen Parsers for Q-Chem input file generation and excited state output parsing

*Workflow published: https://doi.org/10.1063/5.0158061*
*and used https://doi.org/10.1021/acs.jctc.4c00005*
