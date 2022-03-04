function [classification_accuracy_avg, classification_accuracy_stderr, classification_accuracy] = vortex_centre_classifer(No_of_Subject,listen_or_answer)
disp(['initiating...'])    
disp(['classifing task conditions (math vs. story) based on vortex centre positions...'])    

%% load preprocessed fMRI data files
disp(['extracting vortex centre positions...'])    

vortex_filt_nega = [];
vortex_filt_pos = [];
for subject = 1:No_of_Subject
    % load preprocessed fMRI data for each subject
    file_name = ['Preprocessed_bandpass_ScaleSigma_29.35_14.93_sub',num2str(subject),'.mat'];      
    load (file_name);   

    data_allsubject = permute(sigBPass,[2,3,4,1]);
    clearvars sigBPass  

    % find phase map of preprocessed fMRI signal 
        for irow = 1:176
            for icol = 1:251
                temp1 = data_allsubject(irow,icol,:);
                phaseSig(irow,icol,:) = angle(hilbert(temp1(:)));

            end
        end
    clearvars temp1

    % produce phase gradient velocity field baesd on phase map (requires
    % function 'anglesubtract.m')
    Vx_flowmap_norm_phase = [];
    Vy_flowmap_norm_phase = [];
    vPhaseX = zeros(size(phaseSig)) ;
    vPhaseY = zeros(size(phaseSig)) ;

    for iTime = 1:size(phaseSig,3)
        for iX = 1:size(phaseSig,1)
            vPhaseX(iX,2:end-1,iTime) = (anglesubtract(phaseSig(iX,3:end,iTime),phaseSig(iX,1:end-2,iTime)))/2 ;
        end
        for iY = 1:size(phaseSig,2)
            vPhaseY(2:end-1,iY,iTime) = (anglesubtract(phaseSig(3:end,iY,iTime),phaseSig(1:end-2,iY,iTime)))/2 ;
        end
    end
        Vx_flowmap_norm_phase = -vPhaseX./sqrt(vPhaseX.^2+vPhaseY.^2) ; % standardize vector length as 1
        Vy_flowmap_norm_phase = -vPhaseY./sqrt(vPhaseX.^2+vPhaseY.^2) ; % standardize vector length as 1

    clearvars vPhaseX  vPhaseY phaseSig 


    % Detect candidate brain vortices via curl value calcualted from phase vector field: anticlockwise vortex only  
    % clockwise and anticlockwise vortices are detected treated seperately

    for time = 1:size(data_allsubject,3)     
        temp1_vx = Vx_flowmap_norm_phase(:,:,time);
        temp1_vy = Vy_flowmap_norm_phase(:,:,time);
        [x,y] = meshgrid(1:size(temp1_vx,2),1:size(temp1_vx,1));
        [curlz(:,:,time) cav(:,:,time)]  = curl(x,y,temp1_vx,temp1_vy);  % curl value of phase vector field    
    end      
    cav_filt_pos = curlz;
    cav_filt_pos(cav_filt_pos<1) = 0;      % detect only anticlockwise vorties
    cav_filt_pos(isnan(cav_filt_pos)) = 0; % filter out data points with curl value < 1

    % vortex core detection based on filtered curl values (requires function
    % pattDetection_V4.m)

    params.minPattTime = 10; % minimum duration of vortex (10 time steps)
    params.minPattSize = 3;  % minimum size of vortex (3x3 matrix)

    [WCentroids,absoluteTime,instantTotalPower,pattSize,patternIdx] = pattDetection_v4(cav_filt_pos,cav_filt_pos,params,0,'CenterOfMass_Amp');
    clearvars instantTotalPower patternIdx pattSize temp1_vx temp1_vy   

    % find x,y coordinates of each anticlockwise vortex at each time point across subjects 
    for ipatt = 1:size(WCentroids,2)
        temp1_xy = round(WCentroids{ipatt});
        temp1_time = absoluteTime{ipatt};
        for t = 1:size(temp1_xy,1)
           vortex_filt_pos_centroids{ipatt,temp1_time(t),subject} = temp1_xy(t,:);    
        end
    end

    % Detect candidate brain vortices via curl value calcualted from phase vector field: clockwise vortex only  
    % clockwise and anticlockwise vortices are detected treated seperately

    for time = 1:size(data_allsubject,3)      
        temp1_vx = Vx_flowmap_norm_phase(:,:,time);
        temp1_vy = Vy_flowmap_norm_phase(:,:,time);
        [x,y] = meshgrid(1:size(temp1_vx,2),1:size(temp1_vx,1));
        [curlz(:,:,time) cav(:,:,time)]  = curl(x,y,temp1_vx,temp1_vy);  % curl value of phase vector field    
    end      
    cav_filt_pos = curlz;
    cav_filt_pos(cav_filt_pos>-1) = 0;     % detect only clockwise vorties
    cav_filt_pos(isnan(cav_filt_pos)) = 0; % filter out voxels points with curl value > -1

    % vortex core detection based on filtered curl values (requires function
    % pattDetection_V4.m)

    params.minPattTime = 10; % minimum duration of vortex (10 time steps)
    params.minPattSize = 3;  % minimum size of vortex (3x3 matrix)

    [WCentroids,absoluteTime,instantTotalPower,pattSize,patternIdx] = pattDetection_v4(cav_filt_pos,cav_filt_pos,params,0,'CenterOfMass_Amp');
    clearvars instantTotalPower patternIdx pattSize temp1_vx temp1_vy

    % find x,y coordinates of each clockwise vortex at each time point across subjects 
    for ipatt = 1:size(WCentroids,2)
        temp1_xy = round(WCentroids{ipatt});
        temp1_time = absoluteTime{ipatt};
        for t = 1:size(temp1_xy,1)
            vortex_filt_nega_centroids{ipatt,temp1_time(t),subject} = temp1_xy(t,:);    
        end
    end 

end


%% combine both clockwise and anticlockwise vortices in the same map across trials and subjects

  
vortex_template_posi_nega_sparse = [];

for subject = 1:No_of_Subject
    for time = 1:size(data_allsubject,3)
    vortex_template = nan(176,251,2);
    % anticlockwise vortices
    if time <= size(vortex_filt_pos_centroids,2)
        for ipatt = 1:size(vortex_filt_pos_centroids,1);
            temp1 = vortex_filt_pos_centroids{ipatt,time,subject};
            if nansum(temp1(:)) == 0; % skip empty cells
                continue
            else
                row = round(temp1(2)); % find the row coordinate of vortex centre
                col = round(temp1(1)); % find the collumn coordinate of vortex centre            
                % include 3x3 nearby voxels as they will
                % always be included in the vortex
                % mark each included voxel as 1 
                vortex_template(row-1:row+1,col-1:col+1,1) = 1; 

            end
        end
    end
    % anticlockwise vortices
    if time <= size(vortex_filt_nega_centroids,2)
        for ipatt = 1:size(vortex_filt_nega_centroids,1);
            temp1 = vortex_filt_nega_centroids{ipatt,time,subject};
            if nansum(temp1(:)) == 0; % skip empty cells
                continue
            else
                row = round(temp1(2));
                col = round(temp1(1));
                % include 3x3 nearby voxels as they will
                % always be included in the vortex
                % mark each included voxel as -1                 
                vortex_template(row-1:row+1,col-1:col+1,2) = -1 ;

            end
        end
    end
    % calculate the combined distribution of both clockwise and
    % anticlockwise vortices for each subject at each time point
    vortex_template_posi_nega = reshape(nansum(vortex_template,3),[176,251]);
    vortex_template_posi_nega(isnan(vortex_template_posi_nega)) = 0; 
    vortex_template_posi_nega_sparse{time,subject} = vortex_template_posi_nega;
    
 end
end

%% Task specific vortex distribution of each subject
disp(['loading task label...'])    
% load task label         
for subject = 1:No_of_Subject 
    % load task label of each subject
    name = dir([pwd,'/Sample Data/Task Label']) ;            
    file_name2 = [name(subject+2).name];                    
    load (file_name2);     
    fullTime_allsubject{subject} = fullTime;                    
    clearvars sigBPass Vx Vy
end


% define parameters
session_duration = 1; % duration of each task session used for analysis
t_duration = size(data_allsubject,3)-session_duration+1; % last time point to be used for analysis
if listen_or_answer == 1 % type of task, 1 for listening tasks, 2 for answering tasks
    t1 = 6; % starting time of a trial session to be analyzed, 6 for listeing tasks
elseif listen_or_answer == 2 % type of task, 1 for listening tasks, 2 for answering tasks
    t1 = 1; % starting time of a trial session to be analyzed, 1 for answering tasks
end

disp(['extracting task specific vortex centre distributions...'])    
vortex_distribution_math = [];
vortex_distribution_story = [];
vortex_distribution_math_matrix = [];
vortex_distribution_story_matrix = [];
for subject = 1:No_of_Subject
        temp1_sparse_accu_1subject = [];
        % extract vortex distribution for each subject
        for time = 1:t_duration
            temp1_sparse_accu_1subject(:,:,time) = full(vortex_template_posi_nega_sparse{time,subject});
        end

        % extract task label for each subject
        temp1 = fullTime_allsubject{subject};
        temp1_time = temp1(:,3);
        % find time points for math tasks (listening or answering)
        if listen_or_answer == 1 % 1 for listening tasks
            count = find(temp1_time==4); % 4 for math listening tasks
        elseif listen_or_answer == 2 % 2 for answering tasks
            count = find(temp1_time==6); % 6 for math answering tasks
        end
        start_end_time_math = temp1(count,:);

        % time points for story tasks (listening or answering)
        if listen_or_answer == 1 % 1 for listening tasks
            count = find(temp1_time==1); % label-1 for story listening tasks
        elseif listen_or_answer == 2 % 2 for answering tasks
            count = find(temp1_time==3); % label-6 for story answering tasks
        end
        start_end_time_story = temp1(count,:);

        % extract task specifc vortex distribution data: math tasks
        for trial = 1:size(start_end_time_math,1);
            temp2_start = start_end_time_math(trial,1)+ t1;
            temp2_end = temp2_start+session_duration-1;
            if temp2_end > t_duration
                break
            end
            vortex_distribution_math{trial,subject} = temp1_sparse_accu_1subject(:,:,temp2_start:temp2_end);
        end

        % extract task specifc vortex distribution data: story tasks
        for trial = 1:size(start_end_time_story,1);
            temp2_start = start_end_time_story(trial,1)+ t1;
            temp2_end = temp2_start + session_duration-1;
            if temp2_end > t_duration
                break
            end
            vortex_distribution_story{trial,subject} = temp1_sparse_accu_1subject(:,:,temp2_start:temp2_end);
        end
end

%% rearrange vortex centre distribution data by trials  
vortex_distribution_story_matrix = [];
vortex_distribution_math_matrix = [];
empty_trial_vortex_math = [];
empty_trial_vortex_story = [];

disp(['rearranging task-specific vortex centre distributions by trials...'])    

for subject = 1:No_of_Subject    
    
    % math task only
    trialNo = 0;
    for trial = 1:size(vortex_distribution_math,1)
        for subject = 1:No_of_Subject
            temp1_math = vortex_distribution_math{trial,subject}; 
            if nansum(temp1_math(:)) == 0 
            empty_trial_vortex_math{trial,subject} = 1;
            end
            if nansum(temp1_math(:)) ~= 0 
            trialNo = trialNo + 1;  
            vortex_distribution_math_matrix(:,:,:,trialNo) = temp1_math;
            end
        end
    end
    
    % story task only  
    trialNo = 0;
    for trial = 1:size(vortex_distribution_story,1)
        for subject = 1:No_of_Subject
            temp1_story = vortex_distribution_story{trial,subject}; 
            if nansum(temp1_story(:)) == 0 
            empty_trial_vortex_story{trial,subject} = 1;
            end
            if nansum(temp1_story(:)) ~= 0 
            trialNo = trialNo + 1;  
            vortex_distribution_story_matrix(:,:,:,trialNo) = temp1_story;
            end
        end
    end       
end
 
%% classify task conditions (math vs. story) based on vortex centre distributions
disp(['iteratively classifying task conditions via vortex centre positions...'])    
classification_accuracy = [];
% find no. of trials for both math and story tasks
trials_math = size(vortex_distribution_math_matrix,4);
trials_story = size(vortex_distribution_story_matrix,4);
    
% iteratively classify task conditions (story vs. math) with random
% permutations of testing and training trials
for iteration = 1:100 
    % 80% random trials for training
    train_trial_math = randperm(trials_math,round(0.8*trials_math)); 
    train_trial_story = randperm(trials_story,round(0.8*trials_story));
    % 20% remaining trials for testing
    test_trial_story = setdiff([1:trials_story],train_trial_story);
    test_trial_math = setdiff([1:trials_math],train_trial_math);     
    % calculate the mean vortex distribution from training trials 
    % for each task condition, as training template
    % each trial lasts 5 time steps
    temp1_story_train_template = nanmean(nanmean(vortex_distribution_story_matrix(:,:,:,train_trial_story),4),3);
    temp1_math_train_template = nanmean(nanmean(vortex_distribution_math_matrix(:,:,:,train_trial_math),4),3);

    % calculate the vortex distribution of each test trial across 5
    % time steps selected from each session
    temp1_math_test = permute(nanmean(vortex_distribution_math_matrix(:,:,:,test_trial_math),3),[1,2,4,3]);
    temp1_story_test = permute(nanmean(vortex_distribution_story_matrix(:,:,:,test_trial_story),3),[1,2,4,3]);

    temp1_train_math = temp1_math_train_template;
    temp1_train_story = temp1_story_train_template;
    
    % classification of task conditions with only math test trials 
    % for each math test trial, calculate its correlation coefficient with  
    % 2 training templates (story and math)
    for trial = 1:size(test_trial_story,2)
        % remove voxels in both test trials and training templates if 0 
        % and nan values are found in either 
        % this step is to ensure the same data size with no nan/0 in 
        % all variables when calculating correlation coefficients
        temp1_test_math = temp1_math_test(:,:,trial);
        temp1_train_math_filt = temp1_train_math(~isnan(temp1_test_math));
        temp1_train_story_filt = temp1_train_story(~isnan(temp1_test_math));
        temp1_test_math = temp1_test_math(~isnan(temp1_test_math));
        temp1_train_math_filt_01 = temp1_train_math_filt./temp1_train_math_filt;
        temp1_train_story_filt_01 = temp1_train_story_filt./temp1_train_story_filt;
        temp1_test_math_01 = temp1_test_math./temp1_test_math;
        filt_notnan_01 =  temp1_train_math_filt_01.*temp1_train_story_filt_01.*temp1_test_math_01;
        temp1_train_story_filt_filt = temp1_train_story_filt(~isnan(filt_notnan_01));
        temp1_train_math_filt_filt = temp1_train_math_filt(~isnan(filt_notnan_01));
        temp1_test_math_filt = temp1_test_math(~isnan(filt_notnan_01));
        
        if size(temp1_test_math_filt(:)) == 1 % remove test trials with only 1 data points
           continue
        end        
        % calculate correlation coefficient between test trial distribution 
        % and 2 trained templates (story and math)
        R1 = corrcoef(temp1_test_math_filt(:),temp1_train_math_filt_filt(:));
        R2_math_math(trial) = R1(2);

        R1 = corrcoef(temp1_test_math_filt(:),temp1_train_story_filt_filt(:));
        R2_math_story(trial) = R1(2);    
    end
    % compare correlation coefficients between math-math and math-story
    % the higher correlation will selected as the predicted task condition 
    R2_math = [R2_math_math' R2_math_story'];
    % use the total number of story trials here to ensure the same numbers  
    % of trials for both math and story tasks are the same, since story 
    % trials are significantly less than math trials
    max_count_math = nan(1,size(test_trial_story,2)); 
    for trial = 1:size(test_trial_story,2)
        if sum(isnan(R2_math(trial,:)),2) == 0 
            if R2_math(trial,1) ~= R2_math(trial,2)
                [a max_count_math(trial)] = nanmax(R2_math(trial,:));
            elseif R2_math(trial,1) == R2_math(trial,2)
                max_count_math(trial) = 2;
            end
        end
    end
    % prediction result will be compared with the ground truth, math-math
    % (1th column) in this case, wrong predictions (or 2nd column) are marked as 0
    max_count_math_correct = max_count_math;
    max_count_math_correct(max_count_math_correct==2) = 0;


    % testing story trials only
    % for each math test trial, calculate its correlation coefficient with  
    % 2 trained templates (story and math)    
    for trial = 1:size(test_trial_story,2)
        % remove voxels in both test trials and training templates if 0 
        % and nan values are found in either 
        % this step is to ensure the same data size with no nan/0 in 
        % all variables when calculating correlation coefficients
        temp1_test_story = temp1_story_test(:,:,trial);
        temp1_train_math_filt = temp1_train_math(~isnan(temp1_test_story));
        temp1_train_story_filt = temp1_train_story(~isnan(temp1_test_story));    
        temp1_test_story = temp1_test_story(~isnan(temp1_test_story));
        temp1_train_math_filt_01 = temp1_train_math_filt./temp1_train_math_filt;
        temp1_train_story_filt_01 = temp1_train_story_filt./temp1_train_story_filt;
        temp1_test_story_01 = temp1_test_story./temp1_test_story;    
        filt_notnan_01 =  temp1_train_math_filt_01.*temp1_train_story_filt_01.*temp1_test_story_01; 
        temp1_train_story_filt_filt = temp1_train_story_filt(~isnan(filt_notnan_01));
        temp1_train_math_filt_filt = temp1_train_math_filt(~isnan(filt_notnan_01));
        temp1_test_story_filt = temp1_test_story(~isnan(filt_notnan_01));
        
        if size(temp1_test_story_filt(:)) == 1 % remove test trials with only 1 data points
           continue
        end 
        % calculate correlation coefficient between test trial distribution 
        % and 2 trained templates (story and math)
        R1 = corrcoef(temp1_test_story_filt(:),temp1_train_math_filt_filt(:));
        R2_story_math(trial) = R1(2);

        R1 = corrcoef(temp1_test_story_filt(:),temp1_train_story_filt_filt(:));
        R2_story_story(trial) = R1(2);    

    end
    % compare correlation coefficients between math-math and math-story
    % the higher correlation will selected as the predicted task condition     
    R2_story = [ R2_story_story' R2_story_math'];
    max_count_story = nan(1,size(test_trial_story,2));
    for trial = 1:size(test_trial_story,2)
        if sum(isnan(R2_story(trial,:)),2) == 0 
            if R2_story(trial,1) ~= R2_story(trial,2)
                [a max_count_story(trial)] = nanmax(R2_story(trial,:));
            elseif R2_story(trial,1) == R2_story(trial,2)
                max_count_story(trial) = 2;
            end
        end
    end
    % prediction result will be compared with the ground truth, math-math
    % (1th column) in this case, wrong predictions (or 2nd column) are marked as 0    
    max_count_story_correct = max_count_story;
    max_count_story_correct(max_count_story_correct==2) = 0;    

    % combined prediction result of both math and story tasks to calculate
    % prediction accuracy as % of all predictions made, for each iteration 
    max_count_correct = [max_count_story_correct max_count_math_correct];
    max_count_correct_notnan = find(~isnan(max_count_correct));
    max_count_correct_notnan_01 = max_count_correct_notnan./max_count_correct_notnan;
    classification_accuracy(:,iteration) = nansum(max_count_correct(:))./nansum(max_count_correct_notnan_01(:)); % 40 trials before, now depends on how many non-nan results
end
% calculate the mean, stadnard deviation and standard errors of 
% classification accuracy across iterations     
classification_accuracy_avg = nanmean(classification_accuracy,2);
classification_accuracy_std = nanstd(classification_accuracy(:));
classification_accuracy_stderr = classification_accuracy_std./sqrt(iteration);

disp(['completing process...'])    

end
    