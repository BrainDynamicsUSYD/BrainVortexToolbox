# BrainVortexToolbox
MATLAB toolbox to automatically detect and analyse rotational wave patterns in the fMRI data, developed by Dr. Pulin Gong's group at University of Sydney. 

## Instructions for use
System dependencies: N/A <br />
Software version: MATLAB 2016b and above (has been tested on MATLAB 2016b, 2020b, and 2020a) <br />
Hardware requirement: N/A

Data format: fMRI data in standard CIFTI grayordinate space, comprising 32K cortical vertices
Data tested: 100 subjects from the Human Connectomb Project (HCP), [https://db.humanconnectome.org/data/projects/HCP_1200]

### Launch: <br />
Sample data is not avaialble in GitHub due to size limitation, please download from the HCP site link above. No subject ID was given as all subjects should be selected randomly from a corhort of 1200 subjects.

Please download all folders from this repository and allocate the raw fMRI data files (i.e., 'tfMRI_LANGUAGE_LR_Atlas.dtseries.nii') and structural data files (i.e.,'L.flat.32k_fs_LR.surf.gii') downloaded from HCP database into subfolders named by the subject ID under 'Raw Data' and 'Data Pos' folders, repectively. For task state data, please also download the task label files for each subject (witihn the 'EVs' subfolder next to the raw fMRI data file).

For example, the data under subject ID 100206 recorded during a language task should be allocated in the following folders: <br />
Raw fMRI data: '/main_folder/Sample Data/Language Task Original 100 sub/Raw Data/100206/tfMRI_LANGUAGE_LR_Atlas.dtseries.nii'; <br />
Structural data: '/main_folder/Sample Data/Language Task Original 100 sub/Data Pos/100206/L.flat.32k_fs_LR.surf.gii'. <br />
Task label files:  '/main_folder/Sample Data/Language Task Original 100 sub/Raw Data/100206/EVs/present_math.txt'; <br />

Run 'Results_main.m' in matlab to generate key results and figures. 



### Main function: 
Results_main.m (directory to all sub-functions, from preprocessing to spiral detection to sprial-based analysis)

### Subfunctions:

Preprocessing_main.m (pre-preocessing of raw fMRI data with/without spatiotemporal bandpassfilters) <br />
load_fMRI.m  <br />
preproc_fRMI.m  <br />
spaceFreq_fMRI.m  <br />
bandpa_fMRI.m  <br />
ft_read_cifti.m  <br />
gifti.m  <br />

spiral_detection_surfilt.m (spiral detection and filtering based on the statistical test against the null model) <br />
anglesubtract.m <br />
pattDetection_v4.m <br />

spiral_distribution_zscore_map_speed_duration_radius_count.m (spiral distribution z-score map, propagation speed, duration, radius and spiral count) <br />

distance_vs_amplitude.m (relationship between fMRI amplitude and its distance to spiral centres) <br />

spiral_interaction_statistics.m (% proportion of spiral-spiral interactions among three interaction types) <br />

TaskLabel_Extract.m (extract task labels of each subject) <br />

task_specific_spiral_distribution.m (trial-averaged sprial distribution maps across all subjects under each task condition) <br />

spiral_contrast_significance.m (spiral contrast significance maps between different task conditions and their density in 7 functional networks) <br />

spiral_classifer_language_task.m (sprial centre-based linear classifier) <br />

task_evoked_unfiltered_fMRI_signal.m (spiral detection based on task-evoked trial-averaged unfiltered fMRI signals) <br />

PCA.m (principle component analysis of the regional rotational dynamics in the task-evoked trial-averaged unfiltered fMRI signals) <br />

region_of_coordination_ROC.m (identification of the ROC)

PhaseVectorAngle_local_classifier.m (local phase vector field-based linear classifier)

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







