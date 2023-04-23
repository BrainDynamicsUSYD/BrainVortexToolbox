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

task_specific_spiral_distribution.m (trial-averaged spiral distribution maps across all subjects under each task condition) <br />

spiral_contrast_significance.m (spiral contrast significance maps between different task conditions and their density in 7 functional networks) <br />

spiral_classifer_language_task.m (spiral centre-based linear classifier) <br />

task_evoked_unfiltered_fMRI_signal.m (spiral detection based on task-evoked trial-averaged unfiltered fMRI signals) <br />

PCA.m (principle component analysis of the regional rotational dynamics in the task-evoked trial-averaged unfiltered fMRI signals) <br />

region_of_coordination_ROC.m (identification of the ROC)

PhaseVectorAngle_local_classifier.m (local phase vector field-based linear classifier)

## Expect output <br />
Key figures in the paper "Interacting spiral wave patterns underlie complex brain dynamics and are related to cognitive processing" will be generated.

figure 2d. Population-combined z-score maps of brain spiral density in the flattened left cortex  <br />
spiral_distribution_zscore_map_speed_duration_radius_count.m 

figure 2e. fMRI amplitude vs. distance from singularity.<br />
distance_vs_amplitude.m

figure 3e. spiral propagation speed distribution.<br />
spiral_distribution_zscore_map_speed_duration_radius_count.m 

figure 4a. spiral interaction type statistics.<br />
spiral_interaction_statistics.m

figure 5 1st and 2nd columns. task-specific spiral distribution maps.<br />
task_specific_spiral_distribution.m

figure 5 3rd column. spiral contrast significance maps.<br />
spiral_contrast_significance.m

figure 6a. spiral contrast significance density across 7 functional networks.<br />
spiral_contrast_significance.m

figure 6b-c. sprial centre-based linear classifer performances.<br />
spiral_classifer_language_task.m

figure7a-c. sprial detected in task-evoked (trial-averaged) unfiltered fMRI signals.<br />
task_evoked_unfiltered_fMRI_signal.m

figure7d. principle component analysis of the regional rotational dynamics in the task-evoked trial-averaged unfiltered fMRI signals.<br />
PCA.m

figure8b & d. identification of the region of coordination (ROC).<br />
region_of_coordination_ROC.m

figure8c. local phase vector field-based linear classifier.<br />
PhaseVectorAngle_local_classifier.m



## Authors

* **YiBen Xu** - yixu4976@uni.sydney.edu.au
* **Xian Long** - [Xian Long](https://github.com/longxian319)
* **Pulin Gong** - pulin.gong@sydney.edu.au







