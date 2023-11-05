# README

### Installation Requirements:

### Project Info/ System Requirements:
- Python version 3.11.
- Linux or Unix based Operating System required to run the project using the instructions below.
- This project was built using Python 3.11 on macOS Sonoma v14.0.
- This project utilizes the concurrent module for threading from pythons standard library. 
- Development machine: 2.8 GHz Quad-Core Intel 1.7 Processor with 16 GB memory. 

### Files in M3D3 directory:

**Input.gmt.txt**: A list of the 12 known genes associated with Fanconi Anemia with there cooresponding loci.

**STRING 1.txt**: A network of connected genes and their edge measurement.

**results.txt**: A subnetwork of connected fa genes, derived from Module 1 Day 3 project. 

**main.py**: A python script that contains the functionality of creating stage 2 subnetworks, creates an instance of the following _*classes*_ and executes thier main functions:

- **bins.py**: Creates and returns two objects (bins and nfaBins). bins contains a dictionary of # edge connections: genes that contain that number of edge connections from STRING 1.txt. nfaBins, a subset of bins object that strictly contains non FA genes.
- **create_individual_nonfa_subnetwork_tread.py**: Creates and returns an individual stage 2 non FA subnetwork.
- **fagenes.py**: Creates and returns an object containing FA genes, generated from Input.gmt.txt
- **nonfagenes.py**: Creates and returns an object containing non FA genes, generated using STRING 1.txt and faGenes object.
- **stage1_subnetworks.py**: Creates and returns 5000 stage 1 FA subnetworks.

#### The following files are all samples, generated from a successful run of the project.

- **bins.json**
- **nfabins.json**
- **stage1_random_subnetworks.json**
- **stage2_random_subnetworks.json**: Note - The current functionality of this project does not create a new version of this file. To create a new version of this file, please uncomment the create_secondary_subnetwork() call within the main function of main.py
- **stage2_random_subnetworks_v1.json**: A supplementary file, referenced in the writeup (M2D2 HW.pdf). This file is not being used in the project. 


### Setup and Configuration:

- Download and extract zip folder containing the two source text files and main.py, into a known directory on host machine.

<hr>

## Run

1. Open an instance of the terminal.

**Note**: If using a conda environment, please activate the conda environment prior to running the main.py file.

2. Navigate to the directory in which the project was extracted to.

    Example: 
        
        cd ~\Desktop\M2D3\

3. Run the main.py file using python interpreter.

    **Note**: This project requires python version 3.11 or higher, to run the main.py file, please assure you are using the correct interpreter. 

    Example:

        python3 main.py
_***The project will take roughly 2 minutes to execute without the call to create_secondary_subnetwork()._

<hr>

## Results:

Upon successful execution, the json files (minus the stage2_random_subnetworks.json, if left commented) will be replaced and the output of this project will print the means of stage 1 and stage 2 subnetworks, along with the resulting empirical p-value.

Note: If desired, uncomment the call to create_secondary_subnetworks() under the main function within the main.py file, to generate a new stage2_random_subnetworks.json file and recieve a new stage2 edge count mean. 


### References:

- [Permutation concepts 1](https://towardsdatascience.com/how-to-use-permutation-tests-bacc79f45749)
- [Permutation concepts 2](https://www.jwilber.me/permutationtest/)
- [concurrent threading module](https://docs.python.org/3/library/concurrent.futures.html)