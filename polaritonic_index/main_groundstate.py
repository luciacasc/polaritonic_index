import rich_click as click 
from pathlib import Path
from sys import exit
from subprocess import run
import numpy as np

from polaritonic_index import read_file
from polaritonic_index import data_analysis as da
from polaritonic_index import utils





__version__ = '1.1.0'




@click.command()
@click.rich_config(help_config=click.RichHelpConfiguration(use_markdown=True, width=60))
@click.option('-f0', '--filenameGround','fName0', required=True, type = str, help = 'Ground state charge density file.')
@click.option('-r', '--remote', default = None, type= str, help = 'Download file from server. Example: username@remote_server:/path .')
@click.option('-rd', '--radius_type', default="covalent", type=str,
              help="Select the appropriate parameter to characterize the spatial distribution of the charge.\n"
                   "Options:\n"
                   "- 'covalent': Retrieves the covalent radius of the element, which represents the distance between two atoms bonded covalently.\n"
                   "- 'van_der_waals': Retrieves the van der Waals radius of the element, which defines the distance of closest approach between non-bonded atoms.\n"
                   "- 'atomic': Retrieves the atomic radius of the element, which is a measure of the size of an atom.")
@click.option('-v', '--verbosity', default=False,  is_flag=True,  required=None, type=bool, help="Print more information.")

def action( fName0:str,  remote:str, verbosity:bool, radius_type:str ):  
    """
    Processes a ground state charge density file to calculate atomic weights and charges.

    This function takes a ground state charge density file, either locally or downloaded from a remote server, and processes it to:

    1. **Extract atomic information:** Reads atomic coordinates, radii, and initial charges from the file.
    2. **Calculate optimization weights:** Determines optimal weights for each atom using a specified radius type (covalent, van der Waals, or atomic) to minimize the error between the number of valence electrons and the computed electronic charge.
    3. **Calculate localized charges:** Computes the localized ground state charge on each atom using the optimized weights.
    4. **Save results:** Writes the calculated weights and charges to a CSV file named 'weights_charges.csv'.

    Args:
        fName0 (str): Path to the ground state charge density file (`.cube` or `.cub` format).
        remote (str, optional): Remote server location to download the file from.
        verbosity (bool, optional): If True, prints additional information during the process.
        radius_type (str, optional): Type of atomic radius to use for calculations ('covalent', 'van_der_waals', or 'atomic').

    Raises:
        ValueError: If the input file format is incorrect or the file does not exist.
        RuntimeError: If there's an error during file download or processing.
    """
      
    verbosity_flag = verbosity

    ###################################################################################################
    ###############################      Ground state density file      ###############################
    ###################################################################################################

    if fName0:
        fName0 = Path(fName0)
        if remote is None:
            if not fName0.exists():
                print(f"The file {fName0} does not exist")
                exit(1)
            if not fName0.suffix == '.cube' and not fName0.suffix == '.cub':
                print("The file name must have the extension .cube or .cub")
                exit(1)
        else:
            if not remote.endswith('/'):
                remote=f'{remote}/'
            a=run(f'scp {remote}{fName0} ./', text=True)
            if a.returncode == 1:
                print(f"The file {fName0} does not exist")
                exit(1)

        # read the file 
        parse_data = read_file.parse_cube_file((f'{fName0}'))
        atomic_info = read_file.read_coordinates(parse_data, radius_type)  

        atom_type = atomic_info[:, 0]
        atomic_coord = atomic_info[:, 1:4]
        radii = atomic_info[:, 4]
        charge = atomic_info[:, 5]


        density_info, vol = read_file.read_density(parse_data) 
    
        density_coord = density_info[:, :3]
        density_per_point = density_info[:, -1]

        # calculate the weigths 
        weights = da.optimization_radius(atom_type, atomic_coord, charge , density_coord , density_per_point, vol, radii, verbosity_flag)
        
        if verbosity_flag:
           sigma = weights*radii
           utils.print_fancy_table(atom_type, weights, radii, sigma)


        # Calculation of the localized ground state charge on each atom, using optimized weights 
        charge_best_fit, _ = da.charge4atom(atomic_coord, density_coord, density_per_point, vol, radii, weights, verbosity = False)

        data = np.column_stack([weights, charge_best_fit])
        np.savetxt('weights_charges.csv', data , delimiter=',', header='Weights,Computed charge')

        print("Data saved in 'weights_charges.csv'. ")



        
   
    


    
######################################################################################################

if __name__ == '__main__':
    action()
    