import numpy as np
import itertools

from polaritonic_index import utils




def parse_cube_file(filename: str):
    """
    Parses the cube file and extracts the necessary information for grid, coordinates and density.
    You can use this function for ground state density, transition density and difference density files.

    Args:
        filename (str): Name of the cube file.

    Returns:
        tuple: Contains:
            - natom (int): number of atoms.
            - origin (tuple): origin coordinates (x0, y0, z0).
            - grid_size (tuple): number of grid points (NX, NY, NZ).
            - steps (tuple): step sizes (stepx, stepy, stepz).
            - coordinates (list): atoms type, charge, coordinates (x,y,z) in Bohr.
            - density (list): Flattened list of density values.
    """
    

    with open(filename, 'r') as f:
        lines = f.read().splitlines()

    if len(lines) < 6:
        raise ValueError("File does not contain enough lines for a valid cube file.")

    try:
        natom, x0, y0, z0 = lines[2].split()
        natom = int(natom)
        x0, y0, z0 = float(x0), float(y0), float(z0)

        NX, stepx, _, _ = lines[3].split()
        NY, _, stepy, _ = lines[4].split()
        NZ, _, _, stepz = lines[5].split()
        NX, NY, NZ = int(NX), int(NY), int(NZ)
        stepx, stepy, stepz = float(stepx), float(stepy), float(stepz)
    except ValueError:
        raise ValueError("Grid and origin information is malformed.")

    if len(lines) < 6 + natom:
        raise ValueError("File does not contain enough lines for atomic coordinates.")

    try:
        coordinates = lines[6:(6 + natom)]
        coordinates = [i.split() for i in coordinates]

        atom_type = list(map(lambda x: int(x[0]), coordinates))
        charge = list(map(lambda x: float(x[1]), coordinates))
        x = list(map(lambda x: float(x[2]), coordinates))
        y = list(map(lambda x: float(x[3]), coordinates))
        z = list(map(lambda x: float(x[4]), coordinates))
    except ValueError:
        raise ValueError("Atomic coordinates are malformed.")

    try:
        density_lines = lines[(6 + natom):]
        density = list(itertools.chain.from_iterable(line.split() for line in density_lines))
        density = [float(d) for d in density]
    except ValueError:
        raise ValueError("Density values are malformed.")

    return natom, (x0, y0, z0), (NX, NY, NZ), (stepx, stepy, stepz), (atom_type, charge, x, y, z), density






def read_coordinates(parsed_data, radius_type)->np.array:
    """
    Extracts the coordinates matrix from the parsed cube file data.

    Args:
        parsed_data (tuple): Data returned from parse_cube_file.

    Returns:
        np.array: Contains: atom_type, cartesian coordinates(x,y,z) in Bohr, covalent radius 
        in Bohr and charge.
    """
    
    _ , _, _ ,_ , coordinates, _ = parsed_data


    atom_type, charge, x, y, z = coordinates
    #radius = list(map(lambda x: utils.retrieve_atom_radius(x), atom_type))
    radius = list(map(lambda x: utils.get_radius(x,radius_type), atom_type))

    # Check if there are any metal atoms
    metals = {
      "Ru": 44,
      "Rh": 45,
      "Pd": 46,
      "Ag": 47,
      "Os": 76,
      "Ir": 77,
      "Pt": 78,
      "Au": 79,
        }
    
    metal_atomic_numbers = set(metals.values())
    if not any(atom in metal_atomic_numbers for atom in atom_type):
        raise ValueError("Error: No metal atoms found in the coordinates.")


    
    coord_matrix = np.column_stack((atom_type, x, y, z, radius, charge))

    return(coord_matrix)









def read_density(parsed_data):
    """
    Extracts the electron density matrix from the parsed cube file data.

    Args:
        parsed_data (tuple): Data returned from parse_cube_file.

    Returns:
        tuple: Contains:
            - density_matrix (np.array): Matrix with columns (x, y, z, density).
            - vol (float): Volume element.
    """

    _, origin, grid_size, steps, _, density = parsed_data
    x0, y0, z0 = origin
    NX, NY, NZ = grid_size
    stepx, stepy, stepz = steps


    x = np.linspace(x0 + stepx, x0 + stepx * NX, NX)
    y = np.linspace(y0 + stepy, y0 + stepy * NY, NY)
    z = np.linspace(z0 + stepz, z0 + stepz * NZ, NZ)

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    density_matrix = np.column_stack((X.ravel(), Y.ravel(), Z.ravel(), density))
    vol = stepx * stepy * stepz

    return density_matrix, vol