# **User manual: Polaritonic Index**
[![DOI](https://zenodo.org/badge/895511493.svg)](https://doi.org/10.5281/zenodo.14236665)


## **TOC**
* [Theoretical background](#theory)
* [Installation](#installation)
	* [System requirements](#requirements)
	* [Installation instructions](#istructions)
* [How the code works](#workflow)
  * [Full Analysis](#full)
  * [Basic Analysis](#basic)
* [Ground State Charge Density Analysis Tool](#ground)
* [Excited State Charge Density Analysis Tool](#exc)
* [References](#ref)



## **Theoretical background** <a name="theory"></a>

**Polaritonic Index (PI)** is a computational tool designed to analyze and classify excitations in hybrid systems composed of a metal nanoparticle and one or more molecules. By examining transition density file, PI can differentiate between metal, molecular and "hybrid" excitations (polaritonic or charge transfer excitations).

**Excitation Analysis**

PI employs a novel algorithm to analyze the nature of excitations in these hybrid systems. This algorithm calculates a **Polaritonic Index** (PI) for each excitation, which represents  the ratio between the transition density on the molecule(s) and the ones on the metal.

The PI is calculated as follows:

<img src="https://latex.codecogs.com/png.image?\LARGE&space;\dpi{100}\bg{white}PI(\omega_k)=\frac{\Delta\rho_{mol}}{\Delta\rho_{cluster}}=\frac{2^{3/2}\frac{\sum_{I=1}^{Q_{mol}}\sum_{i=1}^{N}e^{-\left[\frac{(\vec{r_i}-\vec{R}I)^2}{2\sigma_i^2}\right]}\left|\Delta\rho(\vec{r_i},\omega_k)\right|}{\sum{j}^{Q_{mol}}Z_j\cdot&space;N_j}}{2^{3/2}\frac{\sum_{J=1}^{Q_{cluster}}\sum_{i=1}^{N}e^{-\left[\frac{(\vec{r_i}-\vec{R}j)^2}{2\sigma_j^2}\right]}\left|\Delta\rho(\vec{r_i},\omega_k)\right|}{\sum{j}^{Q_{cluster}}Z_j\cdot&space;N_j}}" title="PI(\omega_k)=\frac{\Delta\rho_{mol}}{\Delta\rho_{cluster}}=\frac{2^{3/2}\frac{\sum_{I=1}^{Q_{mol}}\sum_{i=1}^{N}e^{-\left[\frac{(\vec{r_i}-\vec{R}I)^2}{2\sigma_i^2}\right]}\left|\Delta\rho(\vec{r_i},\omega_k)\right|}{\sum{j}^{Q_{mol}}Z_j\cdot N_j}}{2^{3/2}\frac{\sum_{J=1}^{Q_{cluster}}\sum_{i=1}^{N}e^{-\left[\frac{(\vec{r_i}-\vec{R}j)^2}{2\sigma_j^2}\right]}\left|\Delta\rho(\vec{r_i},\omega_k)\right|}{\sum{j}^{Q_{cluster}}Z_j\cdot N_j}}" />


where:


- **$\Delta \rho_{mol}$** and **$\Delta \rho_{cluster}$** represent the changes in electron density localized on the molecule and the metal cluster, respectively, for a given excitation at frequency $\omega_k$ (pole of Casida's Equation). 
- The index I runs from 1 to the number of atoms belonging at the molecule ($Q_{mol}$),
instead the index i runs over the grid points (from 1 to n).
- $2^{3/2}$ is the normalization term of the gaussian.
- <img src="https://latex.codecogs.com/png.image?\LARGE&space;\dpi{70}\bg{white}\vec{r}_{i}" title="\vec{r}_{i}" /> indicates a grid point, <img src="https://latex.codecogs.com/png.image?\LARGE&space;\dpi{70}\bg{white}\vec{R}_{I}" title="\vec{R}_{I}" /> is the position of the atom.
- $\sigma_i = w_i \cdot rad $, is the standard deviation, determined by multiplying a weight factor  with a radius. The radius, which can be covalent, van der Waals, or atomic, defines the spatial spread of the charge distribution.

In general, a higher PI value indicates a stronger molecular character, while a lower value suggests a metallic character.



**Interpretation of PI Values:**

- $0 \le PI \le 1 $: Purely metallic excitation.
- $1 < PI < 10$: Hybrid excitation with comparable contributions from the metal and molecule. It could be a charge transfer state or polaritonic state.
- $PI \ge 10$: Purely molecular excitation.



## **How the code works**  <a name="workflow"></a>
The code offers two ways to perform analysis: basic and full.

### **Full Analysis** <a name="full"></a>

**Input Files:** 
Requires two files:
- Ground state electronic density file (*.cube or *.cub format)
- Transition density file of the excited state (*.cube or *.cub format)

**Process:**
The `ground_state_analysis` script analyzes the ground state file. It calculates the contribution of each atom to the total electron density. To estimate atomic contributions to the electron density, a Gaussian function is assigned to each atom. The Gaussian width is optimized by minimizing the difference between theoretical and calculated atomic charges. This width is determined by a weight factor and a radius that characterize the spatial distribution of the electronic charges. The user has the flexibility to parameterize the Gaussian using a variety of radii, including covalent, van der Waals, and atomic.
The final output is a CSV file named 'weights_charges.csv'.  The first column of this file lists the optimized weights, while the second column provides the corresponding best-fit values for the ground state charge density. These best-fit values will be used as normalization factors in subsequent `excited_state_analysis` calculations. The normalization factor is calculated based on the total number of "calculated" valence electrons (Z) and the number of atoms of each type (N).

Following the ground state analysis, the generated 'weights_charges.csv' file can be utilized for subsequent analyses. If the file exists, it is directly employed by the code; otherwise, the user is prompted to provide an alternative file. The polaritonic index is then calculated as described in *Equation 1*. In both ground state and excited state analyses, the Gaussian broadening can be parameterized using covalent, atomic, or van der Waals radii. It is essential for the user to ensure consistency in the choice of radii throughout the analysis.




### **Basic Analysis** <a name="basic"></a>
**Input Files:** 
Requires one file:
- Transition density file of the excited state (*.cube or *.cub format).

**Process:**
In basic analysis mode, the ground state analysis, which involves optimizing weights and calculating charges, is bypassed. Instead, the code directly analyzes the transition density (TD) file. The script attempts to load pre-calculated weights from the 'weights_charges.csv' file. If this file is not found, the user is prompted to either specify an alternative file or provide a uniform weight value between 0.5 and 2. This uniform weight will be assigned to all atoms in the system. The normalization term in the calculation will use the total number of valence electrons (Z) and the number of atoms of each type (N).
The polaritonic index is then computed according to *Equation 1*. 


## Installation <a name="installation"></a>

**WARNING**: Installing via pip from PyPI might not work as expected. You may need to manually install the dependencies (matplotlib, pygenmat, rich-click, and threadpoolctl). However, installation directly from the Git repository functions correctly.
### System requirements <a name="requirements"></a>

The software requirements to run the application are:
- Python $\ge$ 3.9.6 and $\le$ 3.13
- Required libraries: `rich_click`, `numpy`, `pandas`, `matplotlib`, `scikit-learn`.
- Information regarding periodic table elements, such as atomic radius, covalent radius, and van der Waals radius, is provided by the Pymatgen library. 
- Custom modules (`read_file`, `data_analysis`, `utils`).
- `.cube` or `.cub` format file for charge density input.

### Installation instructions <a name="istructions"></a>

For Git users:

```bash
git clone https://github.com/luciacasc/polaritonic_index.git
cd polaritonic_index
pip install -r requirements.txt
pip install .
```
From PyPI:

```bash
pip install polaritonic_index
```


## **Ground State Charge Density Analysis Tool** 

### **Overview**
This tool processes ground state charge density data, extracts atomic information, and calculates weights to minimize the error between the number of valence electrons per atom and the computed electronic charge distribution. The results are saved in a CSV file for further analysis.

---

### **Features**
- Parse `.cube` or `.cub` charge density files.
- Retrieve atomic information, including type, coordinates, and radii.
- Calculate optimization weights and electronic charge distribution.
- Support for local and remote files (download via `scp`).
- Output results to a CSV file (`weights_charges.csv`).

---



### **Usage**
### **Command Syntax**
```bash
ground_state_analysis -f0 <filename> [options]
```

### **Options**
Below are details about all available options, their purpose, and examples of usage:

#### **`-f0, --filenameGround`**  
**Description:**  
Specifies the path to the ground state charge density file. This file contains the electronic density data required for processing.  
 

**Example:**  
```bash
ground_state_analysis -f0 charge_density.cube
```

---

#### **`-r, --remote`**  
**Description:**  
Allows the user to specify a remote server location to download the charge density file using `scp`. This is useful if the file is not stored locally. It requires proper authentication.
The file path must follow the format: `username@server:/path/to/file/`.

**Behavior:**  
- If provided, the script will attempt to download the file using `scp`.
- If not provided, the script assumes the file exists locally.  

**Example:**  
```bash
ground_state_analysis -f0 charge_density.cube -r username@server:/path/to/file
```

---

#### **`-rd, --radius_type`**  
**Description:**  
Specifies the parameter used to characterize the spatial charge distribution of atoms in the system.  

**Options:**  
- `covalent` (default):  
  Uses the covalent radius of the atom, representing the average distance between two bonded atoms.  
- `van_der_waals`:  
  Uses the Van der Waals radius, representing the minimum distance of approach for non-bonded atoms.  
- `atomic`:  
  Uses the atomic radius, representing the overall size of the atom.  

**Example:**  
```bash
ground_state_analysis -f0 <filename> -rd van_der_waals
```

---

#### **`-v, --verbosity`**  
**Description:**  
Enables verbose mode to provide detailed output in the console during the execution of the script.  

**Behavior:**  
- If enabled, the script will display additional information, including intermediate results, weights, and radius calculations in a tabular format.
- If disabled, only essential messages and results will be printed.

**Example:**  
```bash
ground_state_analysis -f0 <filename> -v
```

---

## **Output**
The script produces the following outputs:

1. **CSV File (`weights_charges.csv`):**  
   Contains the calculated weights and electronic charge distribution:
   ```
   #Weights,Computed charge
   0.123,1.456
   0.789,2.345
   ...
   ```

2. **Console Logs (if `-v` is enabled):**  
   Displays intermediate calculations and results in a tabular format.

---

## **Error Handling**
- If the specified file does not exist locally:
  ```plaintext
  The file <filename> does not exist
  ```
- If the file format is incorrect:
  ```plaintext
  The file name must have the extension .cube or .cub
  ```
- If the remote file cannot be downloaded:
  ```plaintext
  The file <filename> does not exist
  ```

---

## **Notes**
- The charge density file must be in `.cube` or `.cub` format. The coordinates provided in the file are assumed to be expressed in Bohr.

---




## **Transition Density Analysis Tool**  <a name="exc"></a>

### **Overview**

This tool analyzes transition density files to classify excitations in systems composed of metallic nanoparticles and melecule(s). 
By leveraging automatic or manual clustering methods, the tool calculates Polaritonic Index (PI) to classify excitations as molecular-like, metal-like, or hybrid. 

This utility provides detailed results, visualizations, and options to control the clustering and analysis process.


Here is an overview of how the code works:

1. **Data Pre-processing**: the code starts with pre-processing the input data. The required format for the input file is ‘cube’.\
The system automatically identifies the number and the type of components (molecules and metallic nanoparticles) and segregates them into subsystems. This process facilitates the classification of excitation types by considering one molecule at a time in conjunction with the cluster. The results of this clustering are visualized using 2D and 3D plots, with different components highlighted in distinct colors (using the -p flag). Alternatively, users can manually clustering, entering the line corresponding to the first atom of each component part of your system.               For example, if your system consists of one cluster of 20 atoms and a molecule of 24 atoms, the manual clustering would be `-m 1 -m 20`.

2. **Classification**: the code calculates the Polatitonic Index using the *Equation 1* formula and classifies the excitations distinguishing the nature bewteen metal, molecular and hybrid.


### **Features**
1. **Transition Density (TD) File Analysis**:
   - Reads atomic coordinates and charge densities from `.cube` files.
   - Supports remote file downloads via `scp`.

2. **Clustering**:
   - **Automatic Clustering**: Uses DBSCAN to identify subsystems within the hybrid system.
   - **Manual Clustering**: Allows users to specify the composition of subsystems manually.
   
3. **Excitation Classification**
Computes Polaritonic Index (PI) to classify excitations.

5. **Visualization**:
Plots the clustered subsystems in 2D and 3D space, providing a visual representation of the molecular and metallic domains.

6. **Verbose Output**:
Offers detailed insights into intermediate calculations, normalization terms, and subsystem charges.


---

### **Usage**
### **Command Syntax**
```bash
excited_state_analysis -f <filenameTD> [options]
```

### **Options**
#### **`-f, --filenameTD`**  
**Description**:  
Specifies the path to the transition density (TD) file. This is the primary input file for analyzing excitations.


**Example**:
```bash
excited_state_analysis -f transition_density.cube
```

---

#### **`-r, --remote`**  
**Description**:  
Downloads the specified file(s) from a remote server using `scp`.

**Requirements**:  
- Optional argument.
- Must follow the format: `username@server:/path/to/file`.

**Example**:
```bash
excited_state_analysis -f transition_density.cube -r username@server:/path/to/files/
```

---

#### **`-rd, --radius_type`**  
**Description**:  
Specifies the type of radius to use for charge density calculations.

**Options**:
- `covalent` (default): Uses the covalent radius for bonded atoms.
- `van_der_waals`: Uses the Van der Waals radius for non-bonded atoms.
- `atomic`: Uses the atomic radius for the overall size of the atom.

**Example**:
```bash
excited_state_analysis -f transition_density.cube -rd van_der_waals
```

---

#### **`-m, --manual`**  
**Description**:  
Enables manual clustering by specifying the first atom of each subsystem. Enter the line corresponding to the first atom of each component part of your system.        For example, if your system consists of one cluster of 20 atoms and a molecule of 24 atoms, the manual clustering would be `-m 1 -m 20`.

**Usage**:  
Provide the starting atom indices for each cluster using `-m`.  

**Example**:
```bash
excited_state_analysis -f transition_density.cube -m 1 -m 20
```

---

#### **`-p, --plot`**  
**Description**:  
Generates 3D visualizations of the clustered subsystems.

**Example**:
```bash
excited_state_analysis -f transition_density.cube -p
```

---

#### **`-v, --verbosity`**  
**Description**:  
Enables verbose mode for detailed output during execution.

**Example**:
```bash
excited_state_analysis -f transition_density.cube -v
```

---

## **Outputs**
1. **Polaritonic Index (PI)**:  
   - A classification of excitations into molecular-like, metal-like, or hybrid categories.

3. **Plots (Optional)**:  
   - 3D representations of clustered subsystems.

4. **Console Logs**:  
   - Verbose mode provides intermediate results, charge distributions, and clustering details.

---

## **Error Handling**
- **File Not Found**: If a specified file does not exist:
  ```plaintext
  The file <filename> does not exist
  ```
- **Invalid File Format**: If a file does not have the `.cube` or `.cub` extension:
  ```plaintext
  The file name must have the extension .cube or .cub
  ```
- **Remote Download Failure**: If a remote file cannot be downloaded:
  ```plaintext
  The file <filename> does not exist
  ```

---

## **Example Workflow**
### **Analyze a Local TD File**
```bash
excited_state_analysis -f transition_density.cube
```

### **Analyze a TD File with Remote Download**
```bash
excited_state_analysis-f transition_density.cube -fd difference_density.cube -r username@server:/path/
```

### **Manual Clustering and Plotting**
```bash
excited_state_analysis -f transition_density.cube -m 1 -m 21 -p
```

---

## **Notes**
- The TD files must be in `.cube` or `.cub` format.
- Remote downloads require appropriate server authentication.
- The charge density file must be in `.cube` or `.cub` format. The coordinates provided in the file are assumed to be expressed in Bohr.
- The code is capable of performing calculations on systems composed of one or more molecules and a single metallic nanoparticle, returning one or more indices.

---

## **References** <a name="ref"></a>

- Cascino L., Corni S., and D’Agostino S., **"Revealing the interplay between hybrid and charge-transfer states in polariton chemistry"**, The Journal of Physical Chemistry C 128.7 (2024): 2917-2927.

---
