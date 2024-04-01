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






### Pymatgen Parsers for Q-Chem input file generation and excited state output parsing

*Workflow published: https://doi.org/10.1063/5.0158061*
*and used https://doi.org/10.1021/acs.jctc.4c00005*
