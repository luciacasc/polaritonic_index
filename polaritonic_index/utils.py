import numpy as np
import os
from rich.console import Console
from rich.table import Table
from rich.style import Style
from pymatgen.core.periodic_table import Element
from pymatgen.analysis.local_env import CovalentRadius




def read_weights_or_ask(number_of_atoms: int) -> np.ndarray:
    """
    Reads weights from 'weights_charges.csv' if the file exists. If the file does not exist, prompts the user
    to either select another file or input a valid floating-point number between 0.5 and 2. This number
    is then applied uniformly to all atoms in the system.

    Args:
        number_of_atoms (int): The number of atoms in the system.

    Returns:
        np.ndarray: A NumPy array containing weights for each atom.
    
    Raises:
        ValueError: If the user's input is not a valid float or is outside the allowed range.
    
    Workflow:
    1. If 'weights_charges.csv' exists, use it.
    2. If 'weights_charges.csv' does not exist, ask the user to:
       - Choose a file from which to read the weights and charges.
       - Or enter a number between 0.1 and 5 that will be uniformly applied to all atoms.
    """

    filename = 'weights_charges.csv'
    
    # Check if the default 'weights_charges.csv' file exists
    if os.path.exists(filename):
        try:
            # Read data from the default CSV file
            data = np.genfromtxt(filename, delimiter=',', skip_header=1)
            print("Data read from 'weights_charges.csv'.")

            # Check if the file has exactly two columns
            if data.ndim != 2 or data.shape[1] != 2:
                print("Error: The file 'weights_charges.csv' must have two columns of data: optimized weights and corresponding ground state charges per atom.")
                return None
            
            # Ensure the number of weights matches the number of atoms
            if len(data) != number_of_atoms:
                print(f"Warning: Number of weights in 'weights_charges.csv' ({len(data)}) does not match the number of atoms ({number_of_atoms}).\n(Hint: maybe there's no header in the file?)")
                return None
            return data
        
        except Exception as e:
            print(f"Error reading the file: {e}")
            return None
    else:
        # Ask the user to choose between reading from a file or entering a number
        print("File 'weights_charges.csv' does not exist.")
        while True:
            choice = input("Choose an option:\n1. Provide another file\n2. Enter a number for uniform weights\nEnter 1 or 2: ")
            if choice == '1':
                # Option 1: User chooses another file
                file_choice = input("Please enter the path to the file: ")
                if os.path.exists(file_choice):
                    try:
                        # Read data from the chosen file
                        data = np.genfromtxt(file_choice, delimiter=',', skip_header=1)
                        print(f"Data read from '{file_choice}'.")
                        
                        # Check if the file has exactly two columns
                        if data.ndim != 2 or data.shape[1] != 2:
                            print(f"Error: The file '{file_choice}' must have two columns of data: optimized weights and corresponding ground state charges per atom.")
                            return None

                        # Ensure the number of weights matches the number of atoms
                        if len(data) != number_of_atoms:
                            print(f"Warning: Number of weights in the file ({len(data)}) does not match the number of atoms ({number_of_atoms}).\n(Hint: maybe there's no header in the file?))  ")
                            return None
                        return data
                    
                    except Exception as e:
                        print(f"Error reading the file: {e}")
                        continue  # Ask the user again for valid input
                else:
                    print("The file does not exist. Please try again.")
            elif choice == '2':
                # Option 2: User enters a number
                while True:
                    try:
                        user_input = float(input("Please enter a number between 0.1 and 5: "))
                        if 0.1 <= user_input <= 5:
                            print(f"You entered: {user_input}")
                            return np.full(number_of_atoms, user_input)
                        else:
                            print("Error: The number must be between 0.1 and 5.")
                    except ValueError:
                        print("Error: Please enter a valid floating-point number.")
            else:
                print("Invalid choice. Please enter '1' or '2'.")





def get_radius(atomic_number, radius_type):
    """
    Retrieves the specified radius in Bohr (covalent, van der Waals, or atomic) for an element based on its atomic number.
    The factor 1.89 is the conversion factor angstrom to bohr
    Parameters:
    -----------
    atomic_number : int
        The atomic number of the element for which the radius is to be retrieved.
        
    radius_type : str
        The type of radius to retrieve. Options are:
        - "covalent": Retrieves the covalent radius of the element.
        - "van_der_waals": Retrieves the van der Waals radius of the element.
        - "atomic": Retrieves the atomic radius of the element.
        
    Returns:
    --------
    float or None
        The radius of the element if available, or None if the radius data is not found.
        
    Raises:
    -------
    ValueError
        If an invalid radius_type is specified or if the requested radius is not available for the element.
    """
        

    try:
        # Get the element based on atomic number
        element = Element.from_Z(atomic_number)
        
        # Retrieve the requested radius
        if radius_type == "covalent":
            # Covalent radius
            radius = CovalentRadius.radius.get(element.symbol, None)
            if radius is None:
                raise ValueError(f"Covalent radius not available for element {element.symbol}.")
        
        elif radius_type == "van_der_waals":
            # Van der Waals radius
            radius = element.van_der_waals_radius
            if radius is None:
                raise ValueError(f"Van der Waals radius not available for element {element.symbol}.")
        
        elif radius_type == "atomic":
            # Atomic radius
            radius = element.atomic_radius
            if radius is None:
                raise ValueError(f"Atomic radius not available for element {element.symbol}.")
        
        else:
            raise ValueError("Invalid radius type. Choose from 'covalent', 'van_der_waals', or 'atomic'.")

        return radius
    
    except Exception as e:
        print(f"Error: {e}")
        return None






def create_exp(x, sigma)-> float:
    """
    Calculate the exponential function for a given value `x` and parameter `sigma`.

    The function computes the value of the exponential expression `exp(-0.5 * x / sigma^2)`.

    Args:
        x (float or np.ndarray): The input value(s) for which the exponential is calculated. 
                                 Can be a single float or a NumPy array of floats.
        sigma (float): The parameter `sigma` which scales the input `x`. Must be a positive float.

    Returns:
        float or np.ndarray: The calculated exponential value(s). If `x` is a single float, 
                             the return is a single float. If `x` is a NumPy array, the return 
                             is a NumPy array of floats with the same shape as `x`.
    """

    return np.exp(-0.5 * x / sigma**2, dtype=np.float64) 





def print_fancy_table(atoms, weights, radii, sigma):
    """
    Prints a fancy table in the terminal with columns: Atom, Weights, Radius, Sigma.
    Different atom types are highlighted with different colors in a generalized way.

    Args:
        atoms (list): List of atom indices.
        weights (np.ndarray): Weights corresponding to each atom.
        radii (np.ndarray): Radii corresponding to each atom.
        sigma (np.ndarray): Sigma values corresponding to each atom.
    """
    console = Console()

    # Create a table
    table = Table()

    # Add columns
    table.add_column("", justify="center", style="cyan", no_wrap=True)
    table.add_column("Atom", justify="center", style="cyan", no_wrap=True)
    table.add_column("Weights", justify="center", style="green")
    table.add_column("Radius", justify="center", style="magenta")
    table.add_column("Sigma", justify="center", style="yellow")

    # List of possible background colors
    color_list = [
        "blue", "green", "red", "yellow", "magenta", "cyan", "white", "bright_blue",
        "bright_green", "bright_red", "bright_yellow", "bright_magenta", "bright_cyan"
    ]

    # Dictionary to store the color assigned to each atom type
    atom_color_map = {}
    
    # Variable to keep track of which color to assign next
    color_index = 0

    # Add rows of data
    for i in range(len(atoms)):
        atom = atoms[i]

        # If this atom hasn't been assigned a color yet, assign one
        if atom not in atom_color_map:
            atom_color_map[atom] = Style(bgcolor=color_list[color_index], bold=True)
            color_index = (color_index + 1) % len(color_list)  # Cycle through colors

        # Get the row style based on the atom's assigned color
        row_style = atom_color_map[atom]

        table.add_row(
            f"{i+1}",
            str(atoms[i]),  # Atom type
            f"{weights[i]:.2f}",  # Weights formatted to 2 decimal places
            f"{radii[i]:.2f}",  # Radii formatted to 2 decimal places
            f"{sigma[i]:.2f}",  # Sigma formatted to 2 decimal places
            style=row_style  # Apply the assigned style
        )

    # Print the table
    console.print(table)
