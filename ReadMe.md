# Novel Techniques for Void Filling in Glacier Elevation Change Data Sets
This reposititory contains the main codebase for the paper 
[_Novel Techniques for Void Filling in Glacier Elevation Change Data Sets_](https://doi.org/10.3390/rs12233917 "https://doi.org/10.3390/rs12233917")

This repository allows to reproduce the main experiments of the paper.

## Contents
* 1. [Installation](#Installation)
* 2. [Dataset](#Dataset)
* 3. [Inpainting](#Inpainting)
* 4. [Statistical Analysis](#StatisticalAnalysis)
* 5. [Known Issues](#KnownIssues)

&nbsp;

## <a name="Installation"></a>1. Installation

### 1.0 Prerequisites
To run the code, a installation of *R* and *python* is needed.
We highly recommend to install python via *anaconda* ([_anaconda.org_](anaconda.org)), as it allows to create an isolated 
environment and simplifies the installation of the needed dependencies.


### 1.1 Download Code

To get started, you can download the code from this repository and change to the respective directory using: 
```bash
git clone https://github.com/tseehaus/inpainting-dhdt
cd inpainting-dhdt
```

### 1.2 R dependencies

To run the code, the following R-dependencies are needed:

* raster
* rgeos
* tools
* Rvision


### 1.3 Python dependencies
The main python-dependencies can be installed using *anaconda*:  
```bash
conda create --file environment.yml
conda activate shearinp
```
This will create an environment called *shearinp* and install the dependencies defined in environment.yml. 
Furthermore, the *pyshearlab* library is required. Installation can be done using pip:  
```bash
pip install https://github.com/stefanloock/pyshearlab/archive/master.zip
```

&nbsp;

## <a name="Dataset"></a>2. Dataset
The dataset can be downloaded at [https://doi.pangaea.de/10.1594/PANGAEA.928371](https://doi.pangaea.de/10.1594/PANGAEA.928371).

After downloading and unzipping the dataset, the folder structure should look as follows:
```
    .
    ├── correlation_setup               # Contains input data and results from the correlation setup (Paper: Section 4.2.2)
    ├── center_setup                    # Contains input data and results from the center setup (Paper: Section 4.2.1)
    ├── juneau-setup                    # Contains input data and results from the juneau setup (Paper: Section 4.2.1)
```



&nbsp;

## <a name="Inpainting"></a>3. Inpainting

To perform the inpainting, we provide 2 scripts:

* perform_inpainting_shearlet.py 
* perform_inpainting_telea_ns.R

### 3.1 Shearlet Inpainting

The script *perform_inpainting_shearlet.py*, performs the shearlet inpainting for a fixed parameter setting (SL5, SL6, or SL7). 
To run it, the following parameters need to be specified:

* input file
* mask file
* ice mask file
* output folder 
* parameter setting
* number of processes 

For a detailed description  of its usage, you can run:
```bash
python perform_inpainting_shearlet.py --help
```

### 3.2 Telea and Navier Stokes Inpainting

The script *perform_inpainting_telea_ns.R*, performs the Telea and Navier Stokes  inpainting. The script will perform both, Navier Stokes and Telea inpainting with the radii 2, 5, 8, 10, 15 and 20.
To run it, the following parameters need to be specified:


* output folder
* input file
* mask file
* ice mask file

For a detailed description of its usage, you can run:
```bash
Rscript perform_inpainting_telea_ns.R --help
```

&nbsp;

## <a name="StatisticalAnalysis"></a>4. Statistical Analysis
For the statistical analysis, we provide two scripts:
* perform_stat_analysis_local.R
* perform_stat_analysis_global.R

### 4.1 Large Voids &ndash; Juneau and Center Setup
The script *perform_stat_analysis_local.R* performs the statistical analysis on the subsets *Juneau* and *Center*, as described in Section 4.2.1 in the paper. It reproduces figure 5 and saves the figure to a png file and the data points to a csv file. To run it, the following parameters need to be specified:
* output folder
* inpaint folder
* shape file folder
* input file
* mask file
* ice mask file

For a more detailed description of its usage, you can run:
```bash
Rscript perform_stat_analysis_local.R --help
```

### 4.2 Large Region &ndash; Correlation Setup
The script *perform_stat_analysis_global.R* performs the statistical analysis on the large region, as described in Section 4.2.2 in the paper. It reproduces figure 6 and saves the figure to a png file and the data points (similar to table 3) to a csv file.
To run it, the following parameters need to be specified:

* output folder
* inpaint folder
* input file
* mask file

For a more detailed description of its usage, you can run:
```bash
Rscript perform_stat_analysis_global.R --help
```
&nbsp;

## <a name="KnownIssues"></a>5. Known Issues

### 5.1 Computation Time of Shearlet Inpainting
The shearlet inpainting algorithm is computationally expensive and the software is implemented for cpu-usage. 
Therefore we recommend to only run it on small images (e.g. 512x512 pixels).  
For large images, it may take up to several days to finish execution.

















