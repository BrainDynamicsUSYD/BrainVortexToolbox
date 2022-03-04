function [vortex_distribution_math_matrix,vortex_distribution_story_matrix, vortex_distribution_math, vortex_distribution_story, vortex_distribution_math_trial_averaged_template,vortex_distribution_story_trial_averaged_template, contrast] = Task_specific_vortex_distribution(No_of_Subject,listen_or_answer)
disp(['initiating process...'])    
%% Task specific Vortex distribution of math and story listening tasks 
disp(['calculating task-specific vortex distribution...'])
%% combine data of both clockwise and anticlockwise rotating vortices into one variable 
disp(['combining clockwise and anticlokwise vortex distributions...'])    

for subject = 1:No_of_Subject;
    % load vortex data     
    filename = ['vortex_bandpass_29.35_14.93_subject_',num2str(subject),'.mat'];
    load(filename);
    

    for time = 1:size(vortex_filt_pos,2)
        vortex_template = nan(176,251,2);

        % combine data of clockwise and anticlockwise rotating
        % vortices into one variable   

        % anticlockwise vortex only
        for ipatt = 1:size(vortex_filt_pos,1); % anticlockwise vortex
            temp1 = vortex_filt_pos{ipatt,time};
            if nansum(temp1(:)) == 0;
               continue
            else
                temp_full = full(temp1);
                [row col] = find(temp_full);
                for i = 1:size(row,1)
                   vortex_template(row(i),col(i),1) = temp_full(row(i),col(i)) ;
                end
            end
        end

        % clockwise vortex only
        if time <= size(vortex_filt_nega,2); % ensure time not exceeding the data length of vortex_filt_nega
                for ipatt = 1:size(vortex_filt_nega,1);
                    temp1 = vortex_filt_nega{ipatt,time};
                    if nansum(temp1(:)) == 0;
                        continue
                    else
                        temp_full = full(temp1);
                        [row col] = find(temp_full);

                        for i = 1:size(row,1)
                           vortex_template(row(i),col(i),2) = temp_full(row(i),col(i)) ;
                        end
                    end
                end
        end
            % calculate the combined vortex distribution by summing both clockwise and anticlockwise
            % vortex distributions
            vortex_template_posi_nega = reshape(nansum(vortex_template,3),[176,251]);
            vortex_template_posi_nega(isnan(vortex_template_posi_nega)) = 0; 
            vortex_template_posi_nega_sparse{time,subject} = vortex_template_posi_nega;
    end
end

%% Task specific vortex distribution of each subject
disp(['loading task labels...'])    
for subject = 1:No_of_Subject 
    % load task label of each subject
    name = dir([pwd,'/Sample Data/Task Label']) ;            
    file_name2 = [name(subject+2).name];                    
    load (file_name2);     
    fullTime_allsubject{subject} = fullTime;                    
    clearvars sigBPass Vx Vy
end
 
% define parameters
session_duration = 5; % duration of each task session used for analysis
t_duration = size(vortex_filt_pos,2)-session_duration+1; % last time point to be used for analysis
if listen_or_answer == 1 % type of task, 1 for listening tasks, 2 for answering tasks
    t1 = 6; % starting time of a trial session to be analyzed, 6 for listeing tasks
elseif listen_or_answer == 2 % type of task, 1 for listening tasks, 2 for answering tasks
    t1 = 1; % starting time of a trial session to be analyzed, 1 for answering tasks
end

disp(['extracting vortex distributions across trials and subjects...'])    
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

        % select task specifc fMRI data: math tasks
        for trial = 1:size(start_end_time_math,1);
            temp2_start = start_end_time_math(trial,1)+ t1;
            temp2_end = temp2_start+session_duration-1;
            if temp2_end > t_duration
                break
            end
            vortex_distribution_math{trial,subject} = temp1_sparse_accu_1subject(:,:,temp2_start:temp2_end);
        end


        % select task specifc fMRI data: story tasks
        for trial = 1:size(start_end_time_story,1);
            temp2_start = start_end_time_story(trial,1)+ t1;
            temp2_end = temp2_start + session_duration-1;
            if temp2_end > t_duration
                break
            end
            vortex_distribution_story{trial,subject} = temp1_sparse_accu_1subject(:,:,temp2_start:temp2_end);
        end
end


%% calculate trial and population-averaged vortex distribution

    vortex_distribution_story_matrix = [];
    vortex_distribution_math_matrix = [];
    empty_trial_vortex_math = [];
    empty_trial_vortex_story = [];

    % rearrange vortex distribution data by trials and combined across
    % subjects
disp(['rearranging vortex distribution data by task trials and combine across subjects...'])    
    
    % math task only
for subject = 1:No_of_Subject    
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
% mean vortex distribution across trials and subjects
    
vortex_distribution_math_trial_averaged_template = nanmean(nanmean(vortex_distribution_math_matrix,4),3);
vortex_distribution_story_trial_averaged_template = nanmean(nanmean(vortex_distribution_story_matrix,4),3); 
    

%% visualization of task-specific trial-averaged vortex distribution 
% of story and math listening/answering tasks, and their contrast
% distribution
disp(['visualization of task-specific vortex distribution...'])    

load('parcellation_template.mat')

% remove vorxels with values of 1 or -1 due to lack of data points
vortex_distribution_story_trial_averaged_template(vortex_distribution_story_trial_averaged_template==1) = nan;
vortex_distribution_math_trial_averaged_template(vortex_distribution_math_trial_averaged_template==1) = nan;
vortex_distribution_story_trial_averaged_template(vortex_distribution_story_trial_averaged_template==-1) = nan;
vortex_distribution_math_trial_averaged_template(vortex_distribution_math_trial_averaged_template==-1) = nan;

% confine the trial averaged vortex distribution within the cortical map
vortex_distribution_math_trial_averaged_template = vortex_distribution_math_trial_averaged_template.*parcellation_template./parcellation_template;
vortex_distribution_story_trial_averaged_template = vortex_distribution_story_trial_averaged_template.*parcellation_template./parcellation_template;

% 2D vortex distribution of math tasks
max = nanmax(vortex_distribution_math_trial_averaged_template(:));
min = nanmin(vortex_distribution_math_trial_averaged_template(:));
% ensure the colour bar centred at 0
if max>min
    max_adj = max;
    min_adj = -1*max;
else
    max_adj = -1*min;
    min_adj = min;
end
figure()
pcolor(vortex_distribution_math_trial_averaged_template)
shading interp
colormap jet
colorbar
caxis([min_adj,max_adj]) % 
if listen_or_answer == 1 % 1 for listening tasks
    title(['tasks-pecific trial-averaged vortex distribution template, math listening'])
elseif listen_or_answer == 2 % 2 for answering tasks
    title(['tasks-pecific trial-averaged vortex distribution template, math answering'])
end

% 2D vortex distribution of story tasks
max = nanmax(vortex_distribution_story_trial_averaged_template(:));
min = nanmin(vortex_distribution_story_trial_averaged_template(:));
% ensure the colour bar centred at 0
if abs(max)>abs(min)
    max_adj = max;
    min_adj = -1*max;
else
    max_adj = -1*min;
    min_adj = min;
end
figure()
pcolor(vortex_distribution_story_trial_averaged_template)
shading interp
colormap jet
colorbar
caxis([min_adj,max_adj])
if listen_or_answer == 1 % 1 for listening tasks
    title(['task-specific trial-averaged vortex distribution template, story listening'])
elseif listen_or_answer == 2 % 2 for answering tasks
    title(['task-specific trial-averaged vortex distribution template, story answering'])
end

% contrast distribution between story and math tasks
contrast = abs(vortex_distribution_story_trial_averaged_template - vortex_distribution_math_trial_averaged_template);
max = nanmax(contrast(:));
min = nanmin(contrast(:));
if max>min
    max_adj = max;
    min_adj = -1*max;
else
    max_adj = -1*min;
    min_adj = min;
end
figure()
pcolor(contrast)
shading interp
colormap jet
colorbar
caxis([0,max_adj])
if listen_or_answer == 1 % 1 for listening tasks
    title(['task-specific vortex distribution contrast, story vs math listening'])
elseif listen_or_answer == 2 % 2 for answering tasks
    title(['task-specific vortex distribution contrast, story vs math answering'])
end


%% task-specific vortex distribution: single voxel vortex count

disp(['visualizing vortex count of a sample voxel...'])    

% define the x,y coordinate of a single sample voxel
x_coordinate = 188;
y_coordinate = 143;

% for each subject, calculate the total vortex count with clockwise and
% anticlockwise rotation directions, as well as number of data points with
% no vortex

% vortex distribution across trials of each subject, math tasks only 
for subject = 1:No_of_Subject
    vortex_distribution_math_matrix = [];
    trialNo = 0;
    for trial = 1:size(vortex_distribution_math,1)
            temp1_math = vortex_distribution_math{trial,subject}; 
            if nansum(temp1_math(:)) == 0 
            empty_trial_vortex_math{trial,subject} = 1;
            end
            if nansum(temp1_math(:)) ~= 0 
            trialNo = trialNo + 1;  
            vortex_distribution_math_matrix(:,:,:,trialNo) = temp1_math;
            end 
    end
    % for each subject, find the vortices of clockwise (-1)and anticlockwise  
    % (1) rotation directions as well as data points with no vortices (0)
    vortex_distribution_math_matrix_1voxel = vortex_distribution_math_matrix(y_coordinate,x_coordinate,:,:);
    temp1 = find(vortex_distribution_math_matrix_1voxel==1); % anticlockwise  votrices (1)
    vortex_distribution_math_matrix_1voxel_poscount(subject) = nansum(temp1./temp1);
    temp2 = find(vortex_distribution_math_matrix_1voxel==-1); % clockwise  votrices (-1)
    vortex_distribution_math_matrix_1voxel_negacount(subject) = nansum(temp2./temp2);
    temp3 = find(vortex_distribution_math_matrix_1voxel==0); % no vortex (0)
    vortex_distribution_math_matrix_1voxel_zerocount(subject) = nansum(temp3./temp3);                
end

% vortex distribution across trials trials of each subject, story tasks only 
for subject = 1:No_of_Subject
        vortex_distribution_story_matrix = [];
        trialNo = 0;
    for trial = 1:size(vortex_distribution_story,1)
            temp1_story = vortex_distribution_story{trial,subject}; 
            if nansum(temp1_story(:)) == 0 
            empty_trial_vortex_story{trial,subject} = 1;
            end
            if nansum(temp1_story(:)) ~= 0 
            trialNo = trialNo + 1;  
            vortex_distribution_story_matrix(:,:,:,trialNo) = temp1_story;
            end
    end   
    % for each subject, find the vortices of clockwise (-1)and anticlockwise  
    % (1) rotation directions as well as data points with no vortices (0)
    vortex_distribution_story_matrix_1voxel = vortex_distribution_story_matrix(y_coordinate,x_coordinate,:,:);
    temp1 = find(vortex_distribution_story_matrix_1voxel==1); % anticlockwise  votrices (1)
    vortex_distribution_story_matrix_1voxel_poscount(subject) = nansum(temp1./temp1);
    temp2 = find(vortex_distribution_story_matrix_1voxel==-1); % clockwise  votrices (-1)
    vortex_distribution_story_matrix_1voxel_negacount(subject) = nansum(temp2./temp2);
    temp3 = find(vortex_distribution_story_matrix_1voxel==0); % no vortex (0)
    vortex_distribution_story_matrix_1voxel_zerocount(subject) = nansum(temp3./temp3);                         
end
           

% calculate the mean and standard error of vortex counts at the sample
% voxel with both rotation directions, across subjects

% story tasks only
vortex_distribution_story_matrix_1voxel_negacount_avg = nanmean(vortex_distribution_story_matrix_1voxel_negacount(:));
vortex_distribution_story_matrix_1voxel_negacount_std = nanstd(vortex_distribution_story_matrix_1voxel_negacount(:));
vortex_distribution_story_matrix_1voxel_negacount_stderr = vortex_distribution_story_matrix_1voxel_negacount_std./sqrt(No_of_Subject);

vortex_distribution_story_matrix_1voxel_poscount_avg = nanmean(vortex_distribution_story_matrix_1voxel_poscount(:));
vortex_distribution_story_matrix_1voxel_poscount_std = nanstd(vortex_distribution_story_matrix_1voxel_poscount(:));
vortex_distribution_story_matrix_1voxel_poscount_stderr = vortex_distribution_story_matrix_1voxel_poscount_std./sqrt(No_of_Subject);

vortex_distribution_story_matrix_1voxel_zerocount_avg = nanmean(vortex_distribution_story_matrix_1voxel_zerocount(:));

vortex_distribution_story_matrix_1voxel_count_avg = [vortex_distribution_story_matrix_1voxel_zerocount_avg vortex_distribution_story_matrix_1voxel_negacount_avg vortex_distribution_story_matrix_1voxel_poscount_avg];

% normalize mean vortex count as % of total number of data poins
vortex_distribution_story_matrix_1voxel_count_avg_norm = vortex_distribution_story_matrix_1voxel_count_avg./nansum(vortex_distribution_story_matrix_1voxel_count_avg(:));

% set the standard error bar
err_low_story = [ vortex_distribution_story_matrix_1voxel_negacount_stderr vortex_distribution_story_matrix_1voxel_poscount_stderr]./nansum(vortex_distribution_story_matrix_1voxel_count_avg(:));
err_high_story = [ vortex_distribution_story_matrix_1voxel_negacount_stderr vortex_distribution_story_matrix_1voxel_poscount_stderr]./nansum(vortex_distribution_story_matrix_1voxel_count_avg(:));

% visualize the distribution with histogram
figure()
b = bar([1:2],vortex_distribution_story_matrix_1voxel_count_avg_norm(2:3));
b.FaceAlpha = 0.5;
hold on
er = errorbar([1:2],vortex_distribution_story_matrix_1voxel_count_avg_norm(2:3),err_low_story,err_high_story);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
if listen_or_answer == 1
    title(['subject-averaged single voxel vortex distribution count, story listening'])
elseif listen_or_answer == 2
    title(['subject-averaged single voxel vortex distribution count, story answering'])
end    
plotName = [{'clockwise'},{'anti-clockwise'}] ;
set(gca,'xticklabel',plotName)

% math tasks only

vortex_distribution_math_matrix_1voxel_negacount_avg = nanmean(vortex_distribution_math_matrix_1voxel_negacount(:));
vortex_distribution_math_matrix_1voxel_negacount_std = nanstd(vortex_distribution_math_matrix_1voxel_negacount(:));
vortex_distribution_math_matrix_1voxel_negacount_stderr = vortex_distribution_math_matrix_1voxel_negacount_std./sqrt(No_of_Subject);

vortex_distribution_math_matrix_1voxel_poscount_avg = nanmean(vortex_distribution_math_matrix_1voxel_poscount(:));
vortex_distribution_math_matrix_1voxel_poscount_std = nanstd(vortex_distribution_math_matrix_1voxel_poscount(:));
vortex_distribution_math_matrix_1voxel_poscount_stderr = vortex_distribution_math_matrix_1voxel_poscount_std./sqrt(No_of_Subject);

vortex_distribution_math_matrix_1voxel_zerocount_avg = nanmean(vortex_distribution_math_matrix_1voxel_zerocount(:));

vortex_distribution_math_matrix_1voxel_count_avg = [vortex_distribution_math_matrix_1voxel_zerocount_avg vortex_distribution_math_matrix_1voxel_negacount_avg vortex_distribution_math_matrix_1voxel_poscount_avg];

% normalize mean vortex count as % of total number of data poins
vortex_distribution_math_matrix_1voxel_count_avg_norm = vortex_distribution_math_matrix_1voxel_count_avg./nansum(vortex_distribution_math_matrix_1voxel_count_avg(:));

% set the standard error bar
err_low_math = [ vortex_distribution_math_matrix_1voxel_negacount_stderr vortex_distribution_math_matrix_1voxel_poscount_stderr]./nansum(vortex_distribution_math_matrix_1voxel_count_avg(:));
err_high_math = [ vortex_distribution_math_matrix_1voxel_negacount_stderr vortex_distribution_math_matrix_1voxel_poscount_stderr]./nansum(vortex_distribution_math_matrix_1voxel_count_avg(:));

% visualize the distribution with histogram
figure()
b = bar([1:2],vortex_distribution_math_matrix_1voxel_count_avg_norm(2:3));
b.FaceAlpha = 0.5;
hold on
er = errorbar([1:2],vortex_distribution_math_matrix_1voxel_count_avg_norm(2:3),err_low_math,err_high_math);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
if listen_or_answer == 1
    title(['subject-averaged single voxel vortex distribution count, math listening'])
elseif listen_or_answer == 2
    title(['subject-averaged single voxel vortex distribution count, math answering'])
end    
plotName = [{'clockwise'},{'anti-clockwise'}] ;
set(gca,'xticklabel',plotName)

disp(['completing process...'])    

end