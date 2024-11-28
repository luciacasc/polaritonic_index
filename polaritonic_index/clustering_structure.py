
def dbscan_model(array, eps=5, min_samples=2):
    """
    Performs DBSCAN (Density-Based Spatial Clustering of Applications with Noise) clustering on a given data array.

    DBSCAN is a clustering algorithm that identifies clusters of high-density data points separated by areas of low density. 
    This function creates and fits a DBSCAN model to the input array, which is assumed to contain coordinate data.

    Args:
        array (numpy.ndarray): An array containing the data points to be clustered. Each row represents an atom, 
            and columns represent features (usually coordinates in this case).
        eps (float, optional): The maximum distance between two points to be considered in a cluster. Defaults to 5.
        min_samples (int, optional): The minimum number of samples in a cluster. Defaults to 2.

    Returns:
        sklearn.cluster.DBSCAN: The fitted DBSCAN clustering model.

    """

    from sklearn.cluster import DBSCAN

    # Create and fit the DBSCAN model with the specified parameters
    clustering_model = DBSCAN(eps=eps, min_samples=min_samples).fit(array)
    return clustering_model


def plot_coordinates(data):
    """
    Generates a scatter plot visualizing the 3D coordinates of each molecule in a provided data dictionary (derived from function "separate_data_by_labels").

    The function extracts all coordinate data from the input dictionary and separates them into x, y, and z components. 
    It then creates a three-subplot figure and plots the coordinates for each molecule in separate subplots, 
    highlighting each molecule with a distinct label.

    Args:
        data (dict): A dictionary containing molecular data. Each key in the dictionary represents a molecule name, 
            and the corresponding value is a list of coordinates for that molecule. Each coordinate is expected to be a list/array 
            containing four elements: [molecule_id, x, y, z].
    """

    import matplotlib.pyplot as plt
    import numpy as np

    # Collecting all coordinates
    all_coordinates = []
    for molecule, coordinates in data.items():
        all_coordinates.extend(coordinates)

    all_coordinates = np.array(all_coordinates)
    x, y, z = all_coordinates[:, 1], all_coordinates[:, 2], all_coordinates[:, 3]

    # Plotting all components together
    plt.figure(figsize=(25, 5))

    plt.subplot(1, 3, 1)
    for molecule, coordinates in data.items():
        coordinates = np.array(coordinates)
        plt.scatter(coordinates[:, 1], coordinates[:, 2], label=f'{molecule} ')
    plt.xlabel('x')
    plt.ylabel('y')
    #plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    plt.subplot(1, 3, 2)
    for molecule, coordinates in data.items():
        coordinates = np.array(coordinates)
        plt.scatter(coordinates[:, 1], coordinates[:, 3], label=f'{molecule} ')
    plt.xlabel('x')
    plt.ylabel('z')
    
    plt.subplot(1, 3, 3)
    for molecule, coordinates in data.items():
        coordinates = np.array(coordinates)
        plt.scatter(coordinates[:, 2], coordinates[:, 3], label=f'{molecule} ')
    plt.xlabel('y')
    plt.ylabel('z')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    plt.suptitle(f'Coordinate Plots for all Components')
    plt.tight_layout(rect=[0, 0, 0.95, 0.95])  # Adjust the layout to accommodate the legend
    plt.show(block = False)
    

def plot_3d_dataset(data:dict):
    """
    Generates a 3D scatter plot visualizing the data points grouped by cluster names.

    This function takes a dictionary containing 3D data points separated into clusters. 
    It creates a 3D scatter plot where each cluster is represented by a distinct color and label.

    Args:
        data (dict): A dictionary where keys are cluster names and values are lists of data points. 
            Each data point is expected to be a list or array with four elements: [cluster_name, X, Y, Z].

    """
    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    for cluster_name, points in data.items():
        
        points = np.array(points)[:, 1:]

        ax.scatter(points[:, 0], points[:, 1], points[:, 2], label=cluster_name)

    ax.set_xlabel('Y')
    ax.set_ylabel('Z')
    ax.set_zlabel('X')  # Invert the axis for a better visualization 
    ax.set_title('3D plot')

    # add a legend
    ax.legend()

    plt.show()


def separate_data_by_labels(data_list, labels):

    """
    Separates data points into clusters based on corresponding labels.

    This function takes a list of data points and a list of corresponding labels, 
    and separates the data points into clusters based on their labels. 
    It performs the following steps:

    1. **Adjusts Labels:** It ensures non-negative labels by subtracting the minimum value from all labels.
    2. **Creates Empty Dictionary:** An empty dictionary `clustered_data` is created to store separated data.
    3. **Iterates Through Data-Label Pairs:** The function iterates through the data points and their corresponding labels.
        - It assigns a cluster name based on the label (`f"molecola_{label}"`).
    4. **Adds Data to Clusters:** The function checks if the cluster name exists in `clustered_data`. 
        - If it exists, the data point is appended to the corresponding cluster list.
        - If it doesn't exist, a new list is created for that cluster name and the data point is added.
    5. **Returns Clustered Data:** Finally, the function returns the dictionary `clustered_data` 
       containing data points separated based on their labels.

    Args:
        data_list (list): A list of data points.
        labels (list): A list of labels corresponding to each data point in `data_list`.

    Returns:
        dict: A dictionary where keys are cluster names (based on labels) and values are lists of data points belonging to that cluster.
        Example: {'metal_cluster':[n_at,x,y,z,r,q], 'molecule_1':[]...., and so on}
    """
    import numpy as np 
    # Find the minimum value of the labels
    min_label = min(labels)
    
    # Add the minimum value to all labels to make them nonnegative
    adjusted_labels = [label - min_label for label in labels]
    
    # Creates an empty dictionary
    clustered_data = {}
    # Iterate through the corresponding data and labels
    for data, label in zip(data_list, adjusted_labels):
        #data_array = np.atleast_2d(data)
        cluster_name = f"molecule_{label}"
        # Add cluster to dictionary
        if cluster_name not in clustered_data:
            clustered_data[cluster_name] = []
        clustered_data[cluster_name].append(data)


    # Check if there are at least 2 structures
    if len(clustered_data) < 2:
        print("Error: Less than 2 structures found. Please check your file or try with manual clusterization (see -m option).")
        return {}
    
    # extract the values from the dict noble_metal
    target_values = set(noble_metal.values())

    new_key_name = "metal_cluster"

    # New dict
    new_clustered_data = {}

    for key, value_list in clustered_data.items():
        concatenated_array = np.concatenate(value_list)
        noble_metal_counts = {metal: np.sum(concatenated_array == metal) for metal in target_values}
        # if statement to check number of metal atoms 
        if any(count >= 5 for count in noble_metal_counts.values()):
            new_clustered_data[new_key_name] = value_list
        else:
            new_clustered_data[key] = value_list

    # Update the dict 
    clustered_data = new_clustered_data

    return clustered_data


noble_metal = {
      "Ru": 44,
      "Rh": 45,
      "Pd": 46,
      "Ag": 47,
      "Os": 76,
      "Ir": 77,
      "Pt": 78,
      "Au": 79,
}


def split_function_manual_mode(molecule:tuple,data,verbosity_flag):
    """
    Splits a data array into clusters based on group start indices provided by -m flag (tuple) and  
    performing metal clustering based on a defined set of noble metals.

    This function takes a tuple, a data array, and an optional verbosity flag. 
    The tuple specifies the starting indices for each group of data points within the data array. 
    The function performs the following steps:

    1. **Splits Data:** It iterates through the molecule tuple and splits the data array accordingly.
        - For the first and last groups, it extracts data from the beginning/end up to the group start index.
        - For intermediate groups, it extracts data between consecutive group start indices.
    2. **Initializes Cluster Dictionary:** An empty dictionary `dict_data` is created to store the split data in clusters.
    3. **Iterates Through Split Data (with Optional Verbosity):**
        - The function iterates through the split data (`split_data`) and assigns a cluster name based on the index (`range(len(split_data))`).
        - If `verbosity_flag` is True, it prints details about the current data cluster and label.
        - It iterates through the elements in the first sub-array of the current data cluster:
            - If `verbosity_flag` is True, it prints the element being checked.
            - It checks if the element is present in the `metalli_nobili` dictionary (containing noble metals).
                - If a noble metal is found, it increments a `metal_count` and prints details (if verbosity is on).
        - Determines the cluster name:
            - If `metal_count` is greater than or equal to 5, the cluster is named "metal_cluster".
            - Otherwise, the cluster name is formatted as "molecola_{label}".
            - If `verbosity_flag` is True, it prints the determined cluster name.
        - Initializes the cluster in `dict_data` if it doesn't exist.
        - Appends all sub-arrays of the current data cluster to the corresponding cluster in `dict_data`.
        - If `verbosity_flag` is True, it prints the updated `dict_data`.
    4. **Returns Clustered Data:** Finally, the function returns the dictionary `dict_data` containing the data points 
       separated into clusters based on group start indices and optionally metal content.

    Args:
        molecule (tuple): A tuple containing group start indices within the data array.
        data (numpy.ndarray): A NumPy array containing the data to be split.
        verbosity_flag (bool, optional): A flag to enable verbosity for debug purposes. 
            If True, the function will print details about the splitting and clustering process. Defaults to False.

    Returns:
        dict: A dictionary where keys are cluster names ("metal_cluster" or "molecola_{label}") and 
              values are lists of sub-arrays (representing data points) belonging to that cluster.
    """

    split_data = []
    for i, group_start in enumerate(molecule):
        if i == 0:
            if len(molecule) > 1:
                next_group_start = molecule[i + 1]
                split_data.append(data[:next_group_start])
            else:
                split_data.append(data[:group_start])
        elif i == len(molecule) - 1:
            split_data.append(data[group_start:])
        else:
            next_group_start = molecule[i + 1]
            split_data.append(data[group_start:next_group_start])
   # Initialize the dictionary to store data clusters
    dict_data = {}

    # Iterate over split_data and their corresponding labels
    for data, label in zip(split_data, range(len(split_data))):
        metal_count = 0  # Reset metal count for each cluster

        # Debug: print the current data substructure and label
        if verbosity_flag==True:
            print('\n Check the clusterization routine: start.')
            print(f"\nProcessing substructure {label}: {data}")

        # Check elements in the first sub-array of the current data cluster
        for element in data[:,0]:
            if element in noble_metal.values():
                metal_count += 1

        if verbosity_flag==True:
            # Debug: print metal count when an element matches
            print(f"Found noble metal: {element}, current metal_count: {metal_count}")

        # Determine the cluster name
        cluster_name = "metal_cluster" if metal_count >= 5 else f"molecule_{label}"
        if verbosity_flag==True:
            # Debug: print the determined cluster name
            print(f"Determined cluster name: {cluster_name}")
            print('Check the clusterization routine: end.')
        # Initialize the cluster in dict_data if not already present
        if cluster_name not in dict_data:
            dict_data[cluster_name] = []
        
        # Append all sub-arrays of the current data cluster to the corresponding cluster in dict_data
        for array_interno in data:
            dict_data[cluster_name].append(array_interno)
        
        

    return dict_data
