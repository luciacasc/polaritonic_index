import numpy as np
import collections
from rich.console import Console
from rich.table import Table
from rich.style import Style

from polaritonic_index import utils 




def optimization_radius(atom_type: np.array, atom_coordinates: np.array, charge : np.array, density_grid : np.array, density_value : np.array, volume: float, radii : np.array, verbosity: bool):
    """
    Optimize the radius by finding the weight that minimizes the difference between theoretical and computed charges.

    Parameters:
    atom_type (np.array): Array of atom types.
    atom_coordinates (np.array): A 2D NumPy array of shape (n_atoms, 3) representing the x, y, z coordinates of the atoms.
    charge (np.array): Array of charges for each atom type (this value is taken from the file).
    density_grid (np.array): A 2D NumPy array of shape (n_points, 3) representing the coordinates of the electron density grid.
    density_value (np.array): A 1D NumPy array of shape (n_points,) containing the electron density values at the grid points.
    volume (float): The volume element used to scale the calculated charge.
    radii (np.array): A 1D NumPy array of shape (n_atoms,) representing the radii of the atoms.
    verbosity (bool): If True, the function will print more info.
  
    Returns:
    new_weights (np.array): The weights that minimize the charge difference.
    """


    weights = np.arange(0.1, 3, 0.1) 
    errors_list = []
    
    
    for i in range(len(weights)):
        charge_calc, _ = charge4atom(atom_coordinates, density_grid, density_value, volume, radii, weights[i], verbosity = False)
        error = abs(charge - charge_calc)
        errors_list.append(error)


    # write table of errors in cvs file
    if verbosity:
        import csv
        csv_file = "errors_vs_weight_table.csv"

        # Write to CSV
        with open(csv_file, mode="w", newline='', encoding="utf-8") as file:
            writer = csv.writer(file)
            # Write header
            header = ["Atom Type"] + [f"Error (w = {weight})" for weight in weights]
            writer.writerow(header)
            
            # Write rows
            for i in range(len(atom_type)):
                row_values = [error[i] for error in errors_list]
                row = [atom_type[i]] + row_values
                writer.writerow(row)

        print(f"Table of errors saved as {csv_file}.")




    # minimum_error =  [min(values) for values in zip(*errors_list)]
    idx_min = [min(enumerate(values), key=lambda x: x[1])[0] for values in zip(*errors_list)]


    new_weights = weights[idx_min]
    
    if verbosity:
        best_charge, _ = charge4atom(atom_coordinates, density_grid, density_value, volume, radii, new_weights, verbosity = True )
        print('\nTotal computed charge: ', round(sum(best_charge),2) , 'vs Total theoretical charge: ', sum(charge))


    return(new_weights)








def charge4atom(atom_coordinates:np.array, density_grid:np.array, density_value:np.array, volume:float, radii: np.array, weights:np.array, verbosity:bool):   
    """
    Calculate the partial charge contribution for a set of atoms based on their coordinates and electron density grid.

    This function computes the partial charge for each atom by evaluating the distance between the atom's coordinates
    and each point in the electron density grid, applying an exponential function to these distances, and then 
    calculating the charge contribution.

    Args:
        atom_coordinates (np.array): A 2D NumPy array of shape (n_atoms, 3) representing the x, y, z coordinates of the atoms.
        density_grid (np.array): A 2D NumPy array of shape (n_points, 3) representing the coordinates of the electron density grid.
        density_value (np.array): A 1D NumPy array of shape (n_points,) containing the electron density values at the grid points.
        volume (float): The volume element used to scale the calculated charge.
        radii (np.array): A 1D NumPy array of shape (n_atoms,) representing the radii of the atoms.
        weight (float): (np.array): A 1D NumPy array of shape (n_atoms,) representing a scaling factor applied to the radii.
        verbosity (bool): If True, the function will print the partial charge contribution.

    Returns:
        np.array: A 1D NumPy array containing the partial charge contributions for each atom.
    """
    x,y,z = atom_coordinates[:, 0], atom_coordinates[:, 1], atom_coordinates[:, 2]
    xg,yg,zg = density_grid[:,0], density_grid[:,1], density_grid[:,2] 
    
    distance = np.array([(xg - x[i])**2 + (yg - y[i])**2 + (zg - z[i])**2 for i in range(len(atom_coordinates))], dtype=np.float64)
    
    sigma = np.array(radii*weights)

    # Initialize an array to store the results
    results = np.zeros_like(distance)
    for i in range(distance.shape[0]):
        results[i] = utils.create_exp(distance[i], sigma[i])


    charge4atom = np.array( 2**(3/2)* (np.dot(results , abs(density_value)))*volume , dtype=np.float64)

    total_charge = sum(charge4atom)
    
    return(charge4atom, total_charge)
    


def normalization(charge: np.array) -> float:
    """
    Calculate the total number of valence electrons in a subsystem based on atomic charges.

    This function takes an array of atomic charges, calculates the number of atoms
    with each unique charge (charge type), and determines the total number of valence 
    electrons for the entire subsystem by multiplying the charge type by the count of atoms
    of that type. The total valence electrons are then returned.

    Parameters:
    ----------
    charge : np.array
        A NumPy array representing the charge values for each atom in the subsystem.

    Returns:
    -------
    float
        The total number of valence electrons for the entire subsystem.
    
    Example:
    --------
    >>> import numpy as np
    >>> charge = np.array([1, 1, 2, 2, 2, 3])
    >>> normalization(charge)
    12.0
    
    Explanation:
    ------------
    In the above example, we have:
    - 2 atoms with a charge of 1 (1 * 2 = 2)
    - 3 atoms with a charge of 2 (2 * 3 = 6)
    - 1 atom with a charge of 3 (3 * 1 = 3)
    The total valence electrons = 2 + 6 + 3 = 11.
    """

    ntype = collections.Counter(charge) 
    valence_electrons_per_atom = [key * value for key, value in ntype.items()]
    total_val_electr_per_subsystem = sum(valence_electrons_per_atom) 

    return (total_val_electr_per_subsystem)





def hybridization_index(dens_metal_array:float,dens_molecule_array:float)->float:
    """This function calculates the hybridization index.

    Args:
        dens_metal_array (float):  Electronic charge of metal cluster.
        dens_molecule_array (float): Electronic charge of each molecule.

    Returns:
        float: Hybridization index 
    """
    
    hi = dens_molecule_array/dens_metal_array
    return(round(hi,2))





def ct_finder(density_grid:np.array, density_value:np.array)->float:
    """
    Calculate the D index, which is a measure of spatial separation of positive and negative densities.

    Args:
        density_grid (np.array): A N x 3 array where each row represents the coordinates (x, y, z) in a 3D space.
        density_value (np.array): density_value (np.array): A 1D NumPy array of shape (n_points,) containing the electron density values at the grid points.

    Returns:
        float: The D index, a scalar value representing the spatial separation between positive and negative density regions.
    """

    rho_p =  np.where(density_value < 0, 0, density_value)
    rho_m = np.where(density_value > 0, 0, density_value)

    tot_p = np.sum(rho_p)
    tot_m = np.sum(rho_m)

    x = density_grid[:,0]
    y = density_grid[:,1]
    z = density_grid[:,2]

    xp= np.dot(x,rho_p)/tot_p
    yp= np.dot(y,rho_p)/tot_p
    zp= np.dot(z,rho_p)/tot_p
    xm= np.dot(x,rho_m)/tot_m
    ym= np.dot(y,rho_m)/tot_m
    zm= np.dot(z,rho_m)/tot_m

    Dx = abs(xp - xm)
    Dy = abs(yp - ym)
    Dz = abs(zp - zm)

    D = (Dx**2 +Dy**2+ Dz**2)**(0.5)

    return(D)









