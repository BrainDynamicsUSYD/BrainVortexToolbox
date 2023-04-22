# BrainVortexToolbox
MATLAB toolbox to automatically detect and analyse rotational wave patterns in the fMRI data, developed by Dr. Pulin Gong's group at University of Sydney. 

## Instructions for use
System dependencies: N/A <br />
Software version: MATLAB 2016b and above (has been tested on MATLAB 2016b, 2020b, and 2021a) <br />
Hardware requirement: N/A

Data format: fMRI data in standard CIFTI grayordinate space, comprising 32K cortical vertices
Data tested: 100 subjects from the Human Connectomb Project (
HCP, https://db.humanconnectome.org/app/template/Login.vm;jsessionid=891FD879A328E1BB3F1B13BAE7655A9E) 

### Launch: <br />
Sample data is not avaialble in GitHub due to size limitation, please download from the HCP site link above.

Create subfolder 'Raw Data' within 'Sample Data' folder, then allocate raw fMRI data into it (language task data file: 'tfMRI_LANGUAGE_LR_Atlas.dtseries.nii', resting data file: 'rfMRI_REST1_LR_Atlas.dtseries.nii';raw fMRI data file should be stored within a subfolder named by the subject ID, I.e., 'Sample Data/Raw Data/100206 (subject ID)/tfMRI_LANGUAGE_LR_Atlas.dtseries.nii').

Create subfolder 'Data Pos' within 'Sample Data' folder, then allocate position file into it (position file: 'L.flat.32k_fs_LR.surf.gii' (left hemisphere), 'R.flat.32k_fs_LR.surf' (right hemisphere); position file should be stored within a subfolder named by the subject ID, I.e., 'Sample Data/Data Pos/100206 (subject ID)/L.flat.32k_fs_LR.surf.gii').

Create subfolder 'Task Label' within 'Sample Data' folder, then allocate language task label file into it (task label file: 'language_label_100206.mat', where 100206 represents subject ID; task label file should be stored directly in the subfolder 'Task Label', I.e.,  'Sample Data/Task Label/language_label_100206.mat').
Run'Main_results.m' in matlab to generate key results and figures. 



### main function: 
Main_results.m (vortex detection)

### subfunctions:

Preprocessing_main.m (pre-preocessing of raw fMRI data with spatiotemporal bandpassfilters) <br />
load_fMRI.m () <br />
preproc_fRMI.m () <br />
spaceFreq_fMRI () <br />
surrogate_fMRI () <br />
bandpa_fMRI () <br />
ft_read_cifti ()  <br />
gifti () <br />

## Expect output <br />
Key figures in the paper "Interacting rotational wave patterns underlie complex brain dynamics and cognitive processing" will be generated.
An example movie (movS1.mov) of detected brain vortices will be generated.

The 'Expected Output Figure' folder contains sample results from two sample subjects.

figure 2. vortex visualization <br />
vortex_detection.m (vortex detection based on curl values of phase vector field) <br />
pattDetection_v4.m (vortex core detection and isolation via cluster analysis) <br />
anglesubtract.m (anglar subtraction) <br />


figure 3. dynamical properties of vortices



figure 4. S2. interaction of brain vortices


figure 5. formation mechanism of brain vortices <br />

vortex_detection_centreonly.m (vortex centre location) <br />
Distance_vs_Amplitude_1subject.m (distance-to-vortex-centre vs. fMRI amplitude, of a single subject) <br />
Distance_vs_Amplitude_Allsubject.m (distance-to-vortex-centre vs. fMRI amplitude, of a all subjects) <br />

figure 6. functional relevance of brain vortices <br />

Task_specific_vortex_distribution.m (task-specific, rotation-specific and population-averaged vortex distribution) <br />
vortex_contrast_distribution_7networks.m (contrast densities between vortex distributions of math and story tasks in 7 functional networks) <br />
parcellation_template.mat (parcellation template of 22 functional areas) <br />
parcellation_template7.mat (parcellation template of 7 functional networks) <br />
vortex_centre_classifer.m (classification accuracy between math and story task conditions based on single-trial moment-by-moment vortex centre distributions)

figure 7. Task-specific activity flow flexibly organized by arrays of brain vortices <br />

Task_specific_phase_vector_field.m (trial-averaged task-specific phase vector field during math listening and answering tasks + 2D heatmap of angle differences between trial averaged phase vector fields of the two tasks + corresponding streamlines) <br />
PhaseVectorAngle_local_classifier.m (2D heatmap of classifiation accuracy between math listening and answering tasks, based on the phase vector field angles of a 3x3 array centred at each target voxel within the ROI) <br />


figure S3. the low dimensional representation of BOLD activities in the vortex


## Authors

* **Xian Long** - [Xian Long](https://github.com/longxian319)
* **YiBen Xu** 
* **Pulin Gong** - pulin.gong@sydney.edu.au







