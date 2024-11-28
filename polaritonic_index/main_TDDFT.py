import rich_click as click 
from pathlib import Path
from sys import exit
from subprocess import run
import numpy as np

from polaritonic_index import utils
from polaritonic_index import read_file
from polaritonic_index import clustering_structure
from polaritonic_index import data_analysis as da



__version__ = '1.1.0'




@click.command()
@click.rich_config(help_config=click.RichHelpConfiguration(use_markdown=True, width=60))
@click.option('-f', '--filenameTD','fName', required=True, type = str, help = 'Transition density file to analyze')
@click.option('-r', '--remote', default = None, type= str, help = 'Download file from server. Example: username@remote_server:/path')
@click.option('-rd', '--radius_type', default="covalent", type=str,
              help="Select the appropriate parameter to characterize the spatial distribution of the charge.\n"
                   "Options:\n"
                   "- 'covalent': Retrieves the covalent radius of the element, which represents the distance between two atoms bonded covalently.\n"
                   "- 'van_der_waals': Retrieves the van der Waals radius of the element, which defines the distance of closest approach between non-bonded atoms.\n"
                   "- 'atomic': Retrieves the atomic radius of the element, which is a measure of the size of an atom.")
@click.option('-m', '--manual', default = [], required=None, multiple=True, type=int, help='Manual clustering involves entering the line corresponding to the first atom of each component part of your system.\
               For example, if your system consists of one cluster of 20 atoms and a molecule of 24 atoms, the manual clustering would be -m 1 -m 20.')
@click.option('-p', '--plot', default=False, is_flag=True, required=None, type=bool, help="Plot the clusterization results to check them")
@click.option('-v', '--verbosity', default=False,  is_flag=True,  required=None, type=bool, help="Print more information")



def action(fName:str, remote:str, radius_type:str, manual, plot, verbosity:bool ):  
    '''    
    Tool for analyzing excitations in hybrid systems using transition density files.

    This script processes cube files to classify excitations in hybrid systems, such as those composed 
    of a metal nanoparticle and one or more molecules. The tool can handle the transition density (TD) 
    files. It performs automatic clustering of atomic data or allows manual input for clustering, and 
    computes polaritonic indicex (PI) to differentiate between 
    molecular-like, metal-like, or hybrid excitations. Optionally, it can also plot the system's 
    clusterization results.

    Usage:
    ------
    The script takes various options for input files, verbosity, and clustering:

    - Transition Density (TD) File: Use the `-f` or `--filenameTD` option to provide the TD file (.cube or .cub).
    - Remote File: Use `-r` or `--remote` to download files from a remote server.
    - Manual Clustering: Use `-m` or `--manual` to manually specify the atomic structure of subsystems.
    - Verbosity: Use `-v` or `--verbosity` for detailed output.
    - Radius type: Use '-rd' or '--radius_type'to select the appropriate parameter to characterize the spatial distribution of the charge.
                   Options:\n
                   - 'covalent': Retrieves the covalent radius of the element, which represents the distance between two atoms bonded covalently.\n
                   - 'van_der_waals': Retrieves the van der Waals radius of the element, which defines the distance of closest approach between non-bonded atoms.\n
                   - 'atomic': Retrieves the atomic radius of the element, which is a measure of the size of an atom.
    - Plotting: Use `-p` or `--plot` to visualize the clustering results.

    Features:
    ---------
    1. **Automatic and Manual Clustering**: 
        - Automatically identifies subsystems using DBSCAN clustering.
        - Optionally allows manual clustering based on user input.

    2. **Excitation Analysis**:
        - Computes the normalized electronic charge for each subsystem.
        - Uses Polaritonic Indix (PI) to classify excitations as:
            - Molecular-like (PI > 10)
            - Metal-like (PI < 1)
            - Hybrid-like (1 <= PI <= 10)

    3. **Plotting**:
        - Visualizes the 3D structure of subsystems after clustering, using either the automatic or manual clustering method.


    Output:
    -------
    - Polaritonic Index to classify the excitation.
    - Clustered subsystem plots, if requested with the `-p` flag.

    
    '''
    
    verbosity_flag = verbosity

    ###################################################################################################
    #######################################         TD file     #######################################
    ###################################################################################################


    fName=Path(fName)
    if remote is None:
        if not fName.exists():
            print(f"The file {fName} does not exist")
            exit(1)
        if not fName.suffix == '.cube' and not fName.suffix == '.cub':
            print("The file name must have the extension .cube or .cub")
            exit(1)
    else:
        if not remote.endswith('/'):
            remote=f'{remote}/'
        a=run(f'scp {remote}{fName} ./', text=True)
        if a.returncode == 1:
            print(f"The file {fName} does not exist")
            exit(1)

    # read the file 
    parse_data = read_file.parse_cube_file((f'{fName}'))
    atomic_info = read_file.read_coordinates(parse_data, radius_type)  

    atomic_coord = atomic_info[:, 1:4]
    radii = atomic_info[:, 4]
    charge = atomic_info[:,-1]


    # set the weights (read from file or ask to user) and normalization term.
    # Normalization term:
        # This 'if' condition determines which normalization term to use for the PI calculation.
        # If the optimization is based on the ground state electron density, the calculated atomic charge is selected; 
        # otherwise, the theoretical charge is chosen.

    weights_charge = utils.read_weights_or_ask(len(atomic_info))

    if weights_charge.ndim !=2:
        weights = weights_charge
        charge_ground = charge     
    else:
        weights = weights_charge[:,0]
        charge_ground = weights_charge[:,1]
        atomic_info[:,-1] = charge_ground # Replace the charge column with atomic information for clustering

        

    if verbosity_flag:
        print('\nNormalization term. If the optimization is based on the ground state electron density, the calculated atomic charge is selected; otherwise, the theoretical charge is chosen.\n')
        print(charge_ground)



    atomic_info = np.column_stack((atomic_info, weights)) # this is necessary to correctly split the weights for each subsystem
    
    
    # clusterization
    if not manual:
        cluster_model = clustering_structure.dbscan_model(atomic_coord) 
        clustered_data = clustering_structure.separate_data_by_labels(atomic_info, cluster_model.labels_) 
        
    else:
        manual_list = clustering_structure.split_function_manual_mode(manual, atomic_info, verbosity)
        


    
    

    # index calculation 
    density_info, vol = read_file.read_density(parse_data) 
    
    density_coord = density_info[:, :3]
    density_per_point = density_info[:, -1]

    dict_charge_system = {}
    if not manual: 
        for key in clustered_data:
            atomic_coord = np.array(clustered_data[key])[:, 1:4]
            radii = np.array(clustered_data[key])[:, 4]
            charge_ground = np.array(clustered_data[key])[:, 5]
            partial_weights = np.array(clustered_data[key])[:,-1]
            
            charge_per_atom, total_charge_subsystem = da.charge4atom(atomic_coord, density_coord, density_per_point, vol, radii, partial_weights, verbosity_flag)
            normalization_term = da.normalization(charge_ground)
            dict_charge_system[key]= total_charge_subsystem/normalization_term
            
            if verbosity_flag:
                print('\nCharge per atom:', key, np.round(charge_per_atom, 5))
                print('\nCharge per subsystem:', key, np.round(total_charge_subsystem, 5))
                print('\nNormalized charge on', key, ':', round(total_charge_subsystem/normalization_term,5) )

    else:
        for key in manual_list:
            atomic_coord = np.array(manual_list[key])[:, 1:4]
            radii = np.array(manual_list[key])[:, 4]
            charge_ground = np.array(manual_list[key])[:,5]
            partial_weights = np.array(manual_list[key])[:,-1]
        
            charge_per_atom, total_charge_subsystem = da.charge4atom(atomic_coord, density_coord, density_per_point, vol, radii, partial_weights, verbosity_flag)
            normalization_term = da.normalization(charge_ground)
            dict_charge_system[key]= total_charge_subsystem/normalization_term
        
            if verbosity_flag:
                print('\nCharge per atom:', key, np.round(charge_per_atom, 5))
                print('\nCharge per subsystem:', key, np.round(total_charge_subsystem, 5))
                print('\nNormalized charge on', key,':', round(total_charge_subsystem/normalization_term,5))


    
    
   
    index_dict = {}
    for key in dict_charge_system:
        if key != 'metal_cluster':
            hybridiz_index = da.hybridization_index(dict_charge_system['metal_cluster'], dict_charge_system[f'{key}'])
            index_dict[key]= hybridiz_index
            
  
    for key, value in index_dict.items():
        if key != 'metal_cluster':
            print(f"Polaritonic Index on \n{key}: {value}.")
            if value > 10:
                print("\nThe excitation should be molecular.")
            elif value < 1:
                print("\nThe excitation should be metal.")
            else:
                print("\nThe excitation should be hybrid. It could be CT or polaritonic state.")
    
 

    
    # plot the system, highlighting the different parts
    if plot:
        if not manual:
            clustering_structure.plot_coordinates(clustered_data)
            clustering_structure.plot_3d_dataset(clustered_data)
           
        else:
            clustering_structure.plot_coordinates(manual_list)
            clustering_structure.plot_3d_dataset(manual_list)



######################################################################################################

if __name__ == '__main__':
    action()
    