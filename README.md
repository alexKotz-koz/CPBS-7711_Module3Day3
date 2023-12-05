# README

### Installation Requirements:

### Project Info/ System Requirements:
- Python version 3.11.
- Linux or Unix based Operating System required to run the project using the instructions below.
- This project was built using Python 3.11 on macOS Sonoma v14.0.
- This project utilizes the concurrent module for threading from pythons standard library. 
- Development machine: 2.8 GHz Quad-Core Intel 1.7 Processor with 16 GB memory. 
- Libraries for visualizaiton:
    - Matplotlib (v3.8.1), NetworkX (v3.2.1), netgraph (v4.12.11)

### Files in M3D3 directory:

**Input.gmt.txt**: A list of the 12 known genes associated with Fanconi Anemia with there cooresponding loci.

**STRING 1.txt**: A network of connected genes and their edge measurement.

**results.txt**: A subnetwork of connected fa genes, derived from Module 1 Day 3 project. 

**M3D3_1000_samples.png**: A visualization of the average gene scores derived from the 1000 randomly generated FA subnetworks (as specified in main.py main() and below). 
- Node: Each node represents a single FA gene
- Node Size: The size of each node corresponds to the average gene score.
- Node Color: The color of each node corresponds to the locus in which the FA gene belongs to.
- Edges: Each edge is derived from the filtered parent network (see discussion above)
- Edges are bundled together to minimize clutter of the image and have no significance.


**main.py**: A python script that contains the functionality of averaging the gene scores and visualizing the average gene scores, creates an instance of the following _*classes*_ and executes thier main functions:

- **fa_utilities.py**: Contains utility functions that can be used for other iterations of this project.
- **py**: Creates and returns an individual stage 2 non FA subnetwork.
- **score_individual_subnet.py**: Creates a gene score for every FA gene, given the individual random FA subnetwork.
- **module2_stage1_subnetworks.py**: Creates and returns 5000 FA subnetworks.

#### The following files are all samples, generated from a successful run of the project.

- **gene.txt**
- **faNetwork.txt**
- **stage1_random_subnetworks.json**


### Setup and Configuration:

- Download and extract zip folder containing the two source text files and main.py, into a known directory on host machine.

<hr>

## Run

1. Open an instance of the terminal.

**Note**: If using a conda environment, please activate the conda environment prior to running the main.py file.

2. Navigate to the directory in which the project was extracted to.

    Example: 
        
        cd ~\Desktop\M3D3\

3. Run the main.py file using python interpreter.

    **Note**: This project requires python version 3.11 or higher, to run the main.py file, please assure you are using the correct interpreter. 

    Example:

        python3 main.py
_***The project will take roughly 52 minutes to execute without adjustments to main.py(), as specified below.***_

<hr>

## Results:

Upon successful execution, a network graph will populate using pythons GUI window that visualizes the average gene scores.

Note: 
- If desired, refer to the comments in the main.py -> main() to run this project on all 5000 random FA subnetworks. 
- **The visualization component takes roughly 5 minutes to execute, depending on the local computer hardware resources**

Google Doc Writeup: https://docs.google.com/document/d/17aR7vwX0a0qanSB3gOwSKhwjZEA5pOgFNOo72rWhN5Y/edit?usp=sharing
