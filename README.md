# BrainVortexToolbox
MATLAB toolbox to automatically detect and analyse rotational wave patterns in the fMRI data, developed by Dr. Pulin Gong's group at University of Sydney. This toolbox includes the code for detecting the rotational wave patterns (brain vortices) in the data, and the code for generating key figures in a manuscript.

## Instructions for use
system dependencies: N/A \n
software version: MATLAB 2016b and above (has been tested on MATLAB 2016b, 2020b, and 2021a) \n
hardware requirement: N/A

data format: fMRI data in standard CIFTI grayordinate space, comprising 32K cortical vertices
data tested: 180 subjects from the Human Connectomb Project (
![HCP](https://db.humanconnectome.org/app/template/Login.vm;jsessionid=891FD879A328E1BB3F1B13BAE7655A9E)) 

Launch: 
1. Add the path of the data (not available due to size limitation, please download from the link above)

2. Type 'main_preprocess' in the matlab command line to generate the preprocessed data (spatial filtering) and detect the vortex statistics.
The preprocessed data and vortex statistics will be saved.

3. Type 'main' in the matlab command line to analyze the vortices and generate key figures.


## main_preprocess.m (vortex detection)
subfunctions:


## main.m (key figures)
figure 2. vortex visualization

figure 3. dynamical properties of vortices

figure 4. S2. interaction of brain vortices

figure 5. formation mechanism of brain vortices

figure S3. the low dimensional representation of BOLD activities in the vortex


## Expect output
Key figures in a manuscript will be generated.
An example movie (movS1.mov) of detected vortices will be gererated.


## Authors

* **Xian Long** - [Xian Long](https://github.com/longxian319)
* **YiBen Xu** 
* **Pulin Gong** - pulin.gong@sydney.edu.au







