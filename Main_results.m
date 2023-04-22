%% Analysis Procedure for main results
clear all
%% parameters 
No_of_Subject = 1; % number of subjects used for analysis
flagSur = 0; %  0 for real data, 1 to generate surrogate data
hemisphere = 1; % 1 for left hemisphere, 2 for right hemisphere
listen_or_answer = 1; % 1 for listening tasks, 2 for answering tasks

addpath(genpath([pwd]))

%% Preprocessing of raw fMRI data from HCP

% first step is to preprocess raw fMRI data using temporal bandpass
% filtering and spatial bandpass filtering (via difference of guassian approach)

    
for subject = 1:No_of_Subject 

    [sigBPass] = Preprocessing_main(subject,hemisphere,flagSur);
    
end

%% Full-size Vortex detection

% Following preprocessing, next is to detect brain vortices based on phase
% vector fields generated from the spatial gradient of phase maps

% We detect the full-sized vortices with an increamentally expanding radius
% and by checking the validity of such vortex using the angle differences between
% the vortex-centre-orgined vector and the phase vector of each voxel

for subject = 1:No_of_Subject
    
    [vortex_filt_nega,vortex_filt_pos,subject] = vortex_detection(subject);
    
end

%% Distance-from-vortex-centre vs. fMRI amplitude
    
% calculate the distance-from-vortex-centre and corresponding amplitude values of each voxel
% within all vortices of each subject

% calculate the vortex-centre-only distributions
for subject = 1:No_of_Subject
    [vortex_filt_nega,vortex_filt_pos,subject] = vortex_detection_centreonly(subject);
end    

%%
% find the relationship between distance-from-vortex-centre and fMRI 
% amplitude of all voxels within all vortices, across subjects
for subject = 1:No_of_Subject    
    [distance_from_centre_amplitude_data_nega_accu,distance_from_centre_amplitude_data_pos_accu] = Distance_vs_Amplitude_1subject(subject);
end      

%%
% after finishing the main calculation, combine the data from each subject
% and plot the mean amplitude (y-axis) corresponding to each unique
% distance-from-vortex-centre values (x-axis) in a line plot

[unique_distanceidx_mm,amplitude_unique_distanceidx_avg_100subjectavg_norm] = Distance_vs_Amplitude_Allsubject(No_of_Subject);

%% Task-specific vortex distribution (trial and population-averaged) 

% task-specific vortex distributions are calculated based on single trial 
% moment-by-moment rotation-direction-specific vortex distributions 
% vortex distribution are displayed in 2D heatmaps for each tasks as well
% as their contrast distribution

% 5 time steps (or 0.72*5 = 3.6s) are selected for both listening and 
% answering tasks for analysis

% histogram of vortex count (clockwise vs. anticlockwise) of single sample
% voxel as % of all data points are also plotted for both math and story
% tasks


[vortex_distribution_math_matrix,vortex_distribution_story_matrix, vortex_distribution_math, vortex_distribution_story, vortex_distribution_math_trial_averaged_template,vortex_distribution_story_trial_averaged_template, contrast] = Task_specific_vortex_distribution(No_of_Subject,listen_or_answer);
% population-averaged vortex contrast density in 7 functional networks
% plot histogram with means and standard errors of vortex contrast density
% distribution (story vs. math) in 7 funtional networks

[contrast_density_parcelavg, contrast_density_parcel_stderr,contrast] = vortex_contrast_distribution_7networks(No_of_Subject,listen_or_answer,vortex_distribution_math,vortex_distribution_story);

%% vortex centre-based linear classifier

% classify task conditions (math vs story) based on single-trial 
% moment-by-moment vortex centre positions
[classification_accuracy_avg, classification_accuracy_stderr, classification_accuracy] = vortex_centre_classifer(No_of_Subject,listen_or_answer);

% due to limited number of trials in sample data, classification performance 
% is significantly reduced 
disp(['iteration-averaged classification performance:'])    
classification_accuracy_avg
disp(['standard error of the averaged classification performance:'])    
classification_accuracy_stderr

%% Task-specific phase vector fields and trial averaged angle difference map between math listening and answering tasks

% calculate trial-averaged phase vector fields for both math
% listening and answering tasks

% calculate and plot the trial-averaged angle differences map between the 
% two task conditions and overlayed with streamlines based on trial
% averaged phase vector fields of either tasks

% plot the angle distribution of all phase vectors at a sample voxel 
% for both math listening and answering tasks with polarhistograms

[vx_math_listen_alltrials,vy_math_listen_alltrials,vx_math_answer_alltrials,vy_math_answer_alltrials] = Task_specific_phase_vector_field(No_of_Subject,listen_or_answer);

%%  phase vector field angle-based local linear classifier

% iteratively (100 iterations with random permutations of test and training 
% trials) classify/predict the task conditions (math listening vs answering) 
% based on single-trial phase vector angles of a local area (3x3 voxels)  
% surrounding each target voxel within a region of interest (ROI)

% mean classification acrruacy across 100 iterations for each target voxel
% within the ROI are ploted in a 2D heatmap to reveal cortical regions with 
% consistent task-specific local phase vector fields at single-trial level

% due to large amount to processing required, this process may take 1+ hour 
% even with sample data 
[classification_accuracy_avg,classification_accuracy_stderr] = PhaseVectorAngle_local_classifier(No_of_Subject,vx_math_listen_alltrials,vy_math_listen_alltrials,vx_math_answer_alltrials,vy_math_answer_alltrials)



