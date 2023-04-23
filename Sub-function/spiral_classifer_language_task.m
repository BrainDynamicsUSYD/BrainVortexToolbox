function [classification_accuracy_avg, classification_accuracy_stderr, classification_accuracy] = spiral_classifer_language_task(flagRest,flagSur,hemisphere,main_folder,No_of_Subject,flagTask)

%% load preprocessed fMRI data files
spiral_template_posi_nega_sparse = [];

for subject = 1:No_of_Subject 
    if flagTask == 1
        foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Spiral Detected'];
        cd(foldername)
        filename = ['Spiral_detected_surfilt_language_task_orig100_LEFT_sub',num2str(subject),'.mat'];
        load(filename) 
    elseif flagTask == 2
        foldername = [main_folder,'/Sample Data/Language Task Additional 100 sub/Spiral Detected'];
        cd(foldername)
        filename = ['spiral_detected_surfilt_language_task_add100_LEFT_sub',num2str(subject),'.mat'];
        load(filename)         
    end
    
% combine both clockwise and anticlockwise spirals in the same map across trials and subjects
    
    x = 1:251;
    y = 1:175;
    [x_grid, y_grid]= meshgrid(x,y);
    radius_threshold = 5;    % mark the 5 radius area surrounding spiral centre with 1 or -1
    for time = 1:size(spiral_filt_pos_real_centreONLY_95perc_extend,2)
        spiral_template_pos = zeros(175,251);
        spiral_template_nega = zeros(175,251);
        if time <= size(spiral_filt_pos_real_centreONLY_95perc_extend,2)
        
    % anticlockwise spirals
        for ipatt = 1:size(spiral_filt_pos_real_centreONLY_95perc_extend,1);
            temp1 = spiral_filt_pos_real_centreONLY_95perc_extend{ipatt,time};
            if nansum(temp1(:)) == 0; % skip empty cells
                continue
            else
            row = temp1(2); % find the row coordinate of spiral centre
            col = temp1(1); % find the row coordinate of spiral centre
            distance_matrix = sqrt((row - y_grid).^2 + (col - x_grid).^2);
            distance_matrix(distance_matrix>radius_threshold) = nan;
            idx_extend = find(~isnan(distance_matrix));
            spiral_template_pos(idx_extend) = 1; % extend the region to 5x5 surrounding the random point and asign the same value as the centre point                
            end
        end
        end
        
    % clockwise spirals
     if time <= size(spiral_filt_nega_real_centreONLY_95perc_extend,2)    
        for ipatt = 1:size(spiral_filt_nega_real_centreONLY_95perc_extend,1);
            temp1 = spiral_filt_nega_real_centreONLY_95perc_extend{ipatt,time};
            if nansum(temp1(:)) == 0; % skip empty cells
                continue
            else
            row = temp1(2); % find the row coordinate of spiral centre
            col = temp1(1); % find the row coordinate of spiral centre
            distance_matrix = sqrt((row - y_grid).^2 + (col - x_grid).^2);
            distance_matrix(distance_matrix>radius_threshold) = nan;
            idx_extend = find(~isnan(distance_matrix));
            spiral_template_nega(idx_extend) = -1; % extend the region to 5x5 surrounding the random point and asign the same value as the centre point      
            end
        end
     end
    % calculate the combined distribution of both clockwise and
    % anticlockwise spirals for each subject at each time point
    spiral_template_posi_nega = spiral_template_nega + spiral_template_pos;
    spiral_template_posi_nega(isnan(spiral_template_posi_nega)) = 0; 
    spiral_template_posi_nega_sparse{time,subject} = spiral_template_posi_nega;
    
 end
end


%% load task label (language task)       

disp(['loading task label...'])    
if flagTask == 1 % language task: original 100 sub
    % load task label of each subject
    foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Task Label'];
    cd(foldername)
    name = dir(pwd) ;   
    file_name2 = ['LanguageOrigTaskLabelAllSubject.mat'];                    
    load (file_name2);     
    for subject = 1:No_of_Subject    
        fullTime_allsubject{subject} = TaskLabel_AllSubject_language_orig{subject};  
    end
elseif flagTask == 2 % language task: additional 100 sub
    % load task label of each subject
    foldername = [main_folder,'/Sample Data/Language Task Additional 100 sub/TaskLabel'];
    cd(foldername)
    name = dir(pwd) ;   
    file_name2 = ['LanguageAddTaskLabelAllSubject.mat'];                    
    load (file_name2);     
    for subject = 1:No_of_Subject    
        fullTime_allsubject{subject} = TaskLabel_AllSubject_language_add{subject};  
    end     
end

%% Task specific spiral distribution of each subject

% define parameters
session_duration = 1; % duration of each task session used for analysis
t_duration = size(spiral_template_posi_nega_sparse,1)-session_duration+1; % last time point to be used for analysis
spiral_distribution_story_listen_matrix = [];
spiral_distribution_story_answer_matrix = []; 
spiral_distribution_math_listen_matrix = [];
spiral_distribution_math_answer_matrix = [];
classification_accuracy = [];

disp(['extracting task specific spiral centre distributions...'])    
spiral_distribution_math = [];
spiral_distribution_story = [];
spiral_distribution_math_matrix = [];
spiral_distribution_story_matrix = [];
for subject = 1:No_of_Subject %1:100
        temp1_sparse_accu_1subject = [];
        % extract spiral distribution for each subject
        for time = 1:t_duration
            temp1 = full(spiral_template_posi_nega_sparse{time,subject});
            if nansum(temp1(:)) ~=0
            temp1_sparse_accu_1subject(:,:,time) = temp1;
            end
        end
        if nansum(temp1_sparse_accu_1subject(:)) == 0
           continue 
        end
        % extract task label for each subject
        temp1 = fullTime_allsubject{subject};
        temp1_time = temp1(:,3);
        % find time points for math tasks (listening or answering)
            count = find(temp1_time==4); % 4 for math listening tasks
            start_end_time_math_listen = temp1(count,:);
            count = find(temp1_time==6); % 6 for math answering tasks
            start_end_time_math_answer = temp1(count,:);

        % time points for story tasks (listening or answering)
            count = find(temp1_time==1); % label-1 for story listening tasks
            start_end_time_story_listen = temp1(count,:);
            count = find(temp1_time==3); % label-6 for story answering tasks
            start_end_time_story_answer = temp1(count,:);

        % extract task specifc spiral distribution data: math tasks
        for trial = 1:size(start_end_time_math_listen,1);
            temp2_start = start_end_time_math_listen(trial,1);
            temp2_end = start_end_time_math_listen(trial,2);;
            temp2_mid = round((temp2_end-temp2_start)./2)+temp2_start;
            temp2_25perc = round((temp2_end-temp2_start).*0.25)+temp2_start;
            temp2_75perc = round((temp2_end-temp2_start).*0.75)+temp2_start;
            t_select = [temp2_start temp2_25perc temp2_mid temp2_75perc temp2_end];
            if temp2_start<1
                continue
            end
            if temp2_end > t_duration
                break
            end
            spiral_distribution_math_listen{trial,subject} = temp1_sparse_accu_1subject(:,:,t_select);
        end
        for trial = 1:size(start_end_time_math_answer,1);
            temp2_start = start_end_time_math_answer(trial,1);                
            temp2_end = start_end_time_math_answer(trial,2);;
            temp2_mid = round((temp2_end-temp2_start)./2)+temp2_start;
            temp2_25perc = round((temp2_end-temp2_start).*0.25)+temp2_start;
            temp2_75perc = round((temp2_end-temp2_start).*0.75)+temp2_start;
            t_select = [temp2_start temp2_25perc temp2_mid temp2_75perc temp2_end];
            if temp2_start<1
                continue
            end
            if temp2_end > t_duration 
                break
            end
            spiral_distribution_math_answer{trial,subject} = temp1_sparse_accu_1subject(:,:,t_select);
        end            
        % extract task specifc spiral distribution data: story tasks
        for trial = 1:size(start_end_time_story_listen,1);
            temp2_start = start_end_time_story_listen(trial,1);                
            temp2_end = start_end_time_story_listen(trial,2);;
            temp2_mid = round((temp2_end-temp2_start)./2)+temp2_start;
            temp2_25perc = round((temp2_end-temp2_start).*0.25)+temp2_start;
            temp2_75perc = round((temp2_end-temp2_start).*0.75)+temp2_start;
            t_select = [temp2_start temp2_25perc temp2_mid temp2_75perc temp2_end];
            if temp2_start<1
                continue
            end                
            if temp2_end > t_duration 
                break
            end
            spiral_distribution_story_listen{trial,subject} = temp1_sparse_accu_1subject(:,:,t_select);
        end
        for trial = 1:size(start_end_time_story_answer,1);
            temp2_start = start_end_time_story_answer(trial,1);                
            temp2_end = start_end_time_story_answer(trial,2);;
            temp2_mid = round((temp2_end-temp2_start)./2)+temp2_start;
            temp2_25perc = round((temp2_end-temp2_start).*0.25)+temp2_start;
            temp2_75perc = round((temp2_end-temp2_start).*0.75)+temp2_start;
            t_select = [temp2_start temp2_25perc temp2_mid temp2_75perc temp2_end];
            if temp2_start<1
                continue
            end                
            if temp2_end > t_duration 
                break
            end
            spiral_distribution_story_answer{trial,subject} = temp1_sparse_accu_1subject(:,:,t_select);
        end            

end

    %% rearrange spiral centre distribution data by trials  

   
    disp(['rearranging task-specific spiral centre distributions by trials...'])    
    for subject = 1:No_of_Subject %1:100    
        % math task only
        trialNo = 0;
        for trial = 1:size(spiral_distribution_math_listen,1)
            for subject = 1:No_of_Subject
                temp1_math = spiral_distribution_math_listen{trial,subject}; 
                if nansum(temp1_math(:)) ~= 0 
                trialNo = trialNo + 1;  
                spiral_distribution_math_listen_matrix(:,:,:,trialNo) = temp1_math;
                end
            end
        end
        trialNo = 0;
        for trial = 1:size(spiral_distribution_math_answer,1)
            for subject = 1:No_of_Subject
                temp1_math = spiral_distribution_math_answer{trial,subject}; 
                if nansum(temp1_math(:)) ~= 0 
                trialNo = trialNo + 1;  
                spiral_distribution_math_answer_matrix(:,:,:,trialNo) = temp1_math;
                end
            end
        end
        % story task only  
        trialNo = 0;
        for trial = 1:size(spiral_distribution_story_listen,1)
            for subject = 1:No_of_Subject
                temp1_story = spiral_distribution_story_listen{trial,subject}; 
                if nansum(temp1_story(:)) ~= 0 
                trialNo = trialNo + 1;  
                spiral_distribution_story_listen_matrix(:,:,:,trialNo) = temp1_story;
                end
            end
        end 
        trialNo = 0;
        for trial = 1:size(spiral_distribution_story_answer,1)
            for subject = 1:No_of_Subject
                temp1_story = spiral_distribution_story_answer{trial,subject}; 
                if nansum(temp1_story(:)) ~= 0 
                trialNo = trialNo + 1;  
                spiral_distribution_story_answer_matrix(:,:,:,trialNo) = temp1_story;
                end
            end
        end
        
    end

    %% classify task conditions (math vs. story) based on spiral centre distributions
    disp(['iteratively classifying task conditions via spiral centre positions...'])    
    % find no. of trials for both math and story tasks
    trials_math_listen = size(spiral_distribution_math_listen_matrix,4);
    trials_math_answer = size(spiral_distribution_math_answer_matrix,4);
    trials_story_listen = size(spiral_distribution_story_listen_matrix,4);
    trials_story_answer = size(spiral_distribution_story_answer_matrix,4);
for t2 = 1:5
    t2
    % iteratively classify task conditions (story vs. math) with random
    % permutations of testing and training trials
    for iteration = 1:100 
        % 80% random trials for training
        train_trial_math_listen = randperm(trials_math_listen,round(0.8*trials_math_listen)); 
        train_trial_math_answer = randperm(trials_math_answer,round(0.8*trials_math_answer)); 
        train_trial_story_listen = randperm(trials_story_listen,round(0.8*trials_story_listen));
        train_trial_story_answer = randperm(trials_story_answer,round(0.8*trials_story_answer));     
        % 20% remaining trials for testing
        test_trial_math_listen = setdiff([1:trials_math_listen],train_trial_math_listen);
        test_trial_math_answer = setdiff([1:trials_math_answer],train_trial_math_answer);
        test_trial_story_listen = setdiff([1:trials_story_listen],train_trial_story_listen);
        test_trial_story_answer = setdiff([1:trials_story_answer],train_trial_story_answer);
        % calculate the mean spiral distribution from training trials 
        % for each task condition, as training template
        % each trial lasts 5 time steps
        temp1_story_listen_train_template = nanmean(nanmean(spiral_distribution_story_listen_matrix(:,:,t2,train_trial_story_listen),4),3);
        temp1_story_answer_train_template = nanmean(nanmean(spiral_distribution_story_answer_matrix(:,:,t2,train_trial_story_answer),4),3);
        temp1_math_listen_train_template = nanmean(nanmean(spiral_distribution_math_listen_matrix(:,:,t2,train_trial_math_listen),4),3);
        temp1_math_answer_train_template = nanmean(nanmean(spiral_distribution_math_answer_matrix(:,:,t2,train_trial_math_answer),4),3);

        % select the spiral distribution of each test trial at 1
        % time step from each session
        temp1_math_listen_test = permute(spiral_distribution_math_listen_matrix(:,:,t2,test_trial_math_listen),[1,2,4,3]);
        temp1_math_answer_test = permute(spiral_distribution_math_answer_matrix(:,:,t2,test_trial_math_answer),[1,2,4,3]);       
        temp1_story_listen_test = permute(nanmean(spiral_distribution_story_listen_matrix(:,:,t2,test_trial_story_listen),3),[1,2,4,3]);
        temp1_story_answer_test = permute(nanmean(spiral_distribution_story_answer_matrix(:,:,t2,test_trial_story_answer),3),[1,2,4,3]);
        
        temp1_train_math_listen = temp1_math_listen_train_template;
        temp1_train_math_answer = temp1_math_answer_train_template;
        temp1_train_story_listen = temp1_story_listen_train_template;
        temp1_train_story_answer = temp1_story_answer_train_template;

        % classification of task conditions with 4 task conditions
        % for each math/story test trial, calculate its correlation coefficient with  
        % 4 training templates (story and math)
        for trial = 1:size(temp1_story_listen_test,3) % use the trial no. of story trials as story trials are significantly less than math trials
            % remove voxels in both test trials and training templates if 0 
            % and nan values are found in either 
            % this step is to ensure the same data size with no nan/0 in 
            % all variables when calculating correlation coefficients
            temp1_test_math_listen = temp1_math_listen_test(:,:,trial);
            temp1_test_math_answer = temp1_math_answer_test(:,:,trial);
            temp1_test_story_listen = temp1_story_listen_test(:,:,trial);
            temp1_test_story_answer = temp1_story_answer_test(:,:,trial);

            % calculate correlation coefficient between test trial distribution 
            % and 4 trained templates (story and math)
            % math listen only
            temp1_test_math_listen_nozeronan_filt = temp1_test_math_listen(temp1_test_math_listen~=0);
            if size(temp1_test_math_listen_nozeronan_filt(:)) == 1 % remove test trials with only 1 data points
               continue
            end                    
            temp1_train_math_listen_nozeronan_filt = temp1_train_math_listen(temp1_test_math_listen~=0);
            temp1_train_math_answer_nozeronan_filt = temp1_train_math_answer(temp1_test_math_listen~=0);
            temp1_train_story_listen_nozeronan_filt = temp1_train_story_listen(temp1_test_math_listen~=0);
            temp1_train_story_answer_nozeronan_filt = temp1_train_story_answer(temp1_test_math_listen~=0);
            
            R1 = corrcoef(temp1_test_math_listen_nozeronan_filt(:),temp1_train_math_listen_nozeronan_filt(:));
            R2_math_listen_math_listen(trial) = R1(2);
            R1 = corrcoef(temp1_test_math_listen_nozeronan_filt(:),temp1_train_math_answer_nozeronan_filt(:));
            R2_math_listen_math_answer(trial) = R1(2);
            R1 = corrcoef(temp1_test_math_listen_nozeronan_filt(:),temp1_train_story_listen_nozeronan_filt(:));
            R2_math_listen_story_listen(trial) = R1(2);            
            R1 = corrcoef(temp1_test_math_listen_nozeronan_filt(:),temp1_train_story_answer_nozeronan_filt(:));
            R2_math_listen_story_answer(trial) = R1(2); 
            
            
            % math answer only
            temp1_test_math_answer_nonanzero_filt = temp1_test_math_answer(temp1_test_math_answer~=0);
            if size(temp1_test_math_answer_nonanzero_filt(:)) == 1 % remove test trials with only 1 data points
               continue
            end                    
            temp1_train_math_listen_nozeronan_filt = temp1_train_math_listen(temp1_test_math_answer~=0);
            temp1_train_math_answer_nozeronan_filt = temp1_train_math_answer(temp1_test_math_answer~=0);
            temp1_train_story_listen_nozeronan_filt = temp1_train_story_listen(temp1_test_math_answer~=0);
            temp1_train_story_answer_nozeronan_filt = temp1_train_story_answer(temp1_test_math_answer~=0);
                        
            R1 = corrcoef(temp1_test_math_answer_nonanzero_filt(:),temp1_train_math_listen_nozeronan_filt(:));
            R2_math_answer_math_listen(trial) = R1(2);
            R1 = corrcoef(temp1_test_math_answer_nonanzero_filt(:),temp1_train_math_answer_nozeronan_filt(:));
            R2_math_answer_math_answer(trial) = R1(2);
            R1 = corrcoef(temp1_test_math_answer_nonanzero_filt(:),temp1_train_story_listen_nozeronan_filt(:));
            R2_math_answer_story_listen(trial) = R1(2);            
            R1 = corrcoef(temp1_test_math_answer_nonanzero_filt(:),temp1_train_story_answer_nozeronan_filt(:));
            R2_math_answer_story_answer(trial) = R1(2);   
            
            % story listen only
            temp1_test_story_listen_nonanzero_filt = temp1_test_story_listen(temp1_test_story_listen~=0);
            if size(temp1_test_story_listen_nonanzero_filt(:)) == 1 % remove test trials with only 1 data points
               continue
            end                    
            temp1_train_math_listen_nozeronan_filt = temp1_train_math_listen(temp1_test_story_listen~=0);
            temp1_train_math_answer_nozeronan_filt = temp1_train_math_answer(temp1_test_story_listen~=0);
            temp1_train_story_listen_nozeronan_filt = temp1_train_story_listen(temp1_test_story_listen~=0);
            temp1_train_story_answer_nozeronan_filt = temp1_train_story_answer(temp1_test_story_listen~=0);
    
            R1 = corrcoef(temp1_test_story_listen_nonanzero_filt(:),temp1_train_math_listen_nozeronan_filt(:));
            R2_story_listen_math_listen(trial) = R1(2);
            R1 = corrcoef(temp1_test_story_listen_nonanzero_filt(:),temp1_train_math_answer_nozeronan_filt(:));
            R2_story_listen_math_answer(trial) = R1(2);
            R1 = corrcoef(temp1_test_story_listen_nonanzero_filt(:),temp1_train_story_listen_nozeronan_filt(:));
            R2_story_listen_story_listen(trial) = R1(2);            
            R1 = corrcoef(temp1_test_story_listen_nonanzero_filt(:),temp1_train_story_answer_nozeronan_filt(:));
            R2_story_listen_story_answer(trial) = R1(2);   
            
            % story answer only
            temp1_test_story_answer_nonanzero_filt = temp1_test_story_answer(temp1_test_story_answer~=0);
            if size(temp1_test_story_answer_nonanzero_filt(:)) == 1 % remove test trials with only 1 data points
               continue
            end                    
            temp1_train_math_listen_nozeronan_filt = temp1_train_math_listen(temp1_test_story_answer~=0);
            temp1_train_math_answer_nozeronan_filt = temp1_train_math_answer(temp1_test_story_answer~=0);
            temp1_train_story_listen_nozeronan_filt = temp1_train_story_listen(temp1_test_story_answer~=0);
            temp1_train_story_answer_nozeronan_filt = temp1_train_story_answer(temp1_test_story_answer~=0);
                
            R1 = corrcoef(temp1_test_story_answer_nonanzero_filt(:),temp1_train_math_listen_nozeronan_filt(:));
            R2_story_answer_math_listen(trial) = R1(2);
            R1 = corrcoef(temp1_test_story_answer_nonanzero_filt(:),temp1_train_math_answer_nozeronan_filt(:));
            R2_story_answer_math_answer(trial) = R1(2);
            R1 = corrcoef(temp1_test_story_answer_nonanzero_filt(:),temp1_train_story_listen_nozeronan_filt(:));
            R2_story_answer_story_listen(trial) = R1(2);            
            R1 = corrcoef(temp1_test_story_answer_nonanzero_filt(:),temp1_train_story_answer_nozeronan_filt(:));
            R2_story_answer_story_answer(trial) = R1(2);              
        end
        
        % compare correlation coefficients between math-math and math-story
        % the higher correlation will selected as the predicted task condition 
        R2_math_listen = [R2_math_listen_math_listen' R2_math_listen_math_answer' R2_math_listen_story_listen' R2_math_listen_story_answer'];
        R2_math_answer = [R2_math_answer_math_listen' R2_math_answer_math_answer' R2_math_answer_story_listen' R2_math_answer_story_answer'];
        R2_story_listen = [R2_story_listen_math_listen' R2_story_listen_math_answer' R2_story_listen_story_listen' R2_story_listen_story_answer'];
        R2_story_answer = [R2_story_answer_math_listen' R2_story_answer_math_answer' R2_story_answer_story_listen' R2_story_answer_story_answer'];
        
        % use the total number of story trials here to ensure the same numbers  
        % of trials for both math and story tasks are the same, since story 
        % trials are significantly less than math trials
%         max_count_math_listen = nan(1,size(test_trial_story_listen,2)); 
        for trial = 1:size(R2_math_listen,1)
            if sum(isnan(R2_math_listen(trial,:)),2) == 0 
                if R2_math_listen(trial,1) ~= R2_math_listen(trial,2)
                    [a max_count_math_listen(trial)] = nanmax(R2_math_listen(trial,:));
                elseif R2_math_listen(trial,1) == R2_math_listen(trial,2)
                    max_count_math_listen(trial) = 2;
                end
            end
            if sum(isnan(R2_math_answer(trial,:)),2) == 0 
                if R2_math_answer(trial,1) ~= R2_math_answer(trial,2)
                    [a max_count_math_answer(trial)] = nanmax(R2_math_answer(trial,:));
                elseif R2_math_answer(trial,1) == R2_math_answer(trial,2)
                    max_count_math_answer(trial) = 1;
                end
            end            
            if sum(isnan(R2_story_listen(trial,:)),2) == 0 
                if R2_story_listen(trial,3) ~= R2_story_listen(trial,4)
                    [a max_count_story_listen(trial)] = nanmax(R2_story_listen(trial,:));
                elseif R2_story_listen(trial,3) == R2_story_listen(trial,4)
                    max_count_story_listen(trial) = 4;
                end
            end  
            if sum(isnan(R2_story_answer(trial,:)),2) == 0 
                if R2_story_answer(trial,3) ~= R2_story_answer(trial,4)
                    [a max_count_story_answer(trial)] = nanmax(R2_story_answer(trial,:));
                elseif R2_story_answer(trial,3) == R2_story_answer(trial,4)
                    max_count_story_answer(trial) = 3;
                end
            end              
        end
        % prediction result will be compared with the ground truth, math-math
        % (1th column) in this case, wrong predictions (or 2nd column) are marked as 0
        % 4 conditions
        max_count_math_listen_correct = max_count_math_listen;
        max_count_math_listen_correct(max_count_math_listen_correct~=1) = 0;        
        max_count_math_answer_correct = max_count_math_answer;
        max_count_math_answer_correct(max_count_math_answer_correct~=2) = 0;
        max_count_math_answer_correct(max_count_math_answer_correct==2) = 1;
        max_count_story_listen_correct = max_count_story_listen;
        max_count_story_listen_correct(max_count_story_listen_correct~=3) = 0;  
        max_count_story_listen_correct(max_count_story_listen_correct==3) = 1;        
        max_count_story_answer_correct = max_count_story_answer;
        max_count_story_answer_correct(max_count_story_answer_correct~=4) = 0; 
        max_count_story_answer_correct(max_count_story_answer_correct==4) = 1; 
        
        % combined prediction result of both math and story tasks to calculate
        % prediction accuracy as % of all predictions made, for each iteration 
        max_count_correct = [max_count_math_listen_correct max_count_math_answer_correct max_count_story_listen_correct max_count_story_answer_correct];
        max_count_correct_notnan = find(~isnan(max_count_correct));
        max_count_correct_notnan_01 = max_count_correct_notnan./max_count_correct_notnan;
        classification_accuracy(iteration,t2) = nansum(max_count_correct(:))./nansum(max_count_correct_notnan_01(:)); % 40 trials before, now depends on how many non-nan results
    end
    % calculate the mean, stadnard deviation and standard errors of 
    % classification accuracy across iterations     
    classification_accuracy_avg(t2) = nanmean(classification_accuracy(:,t2),1);
    classification_accuracy_std(t2) = nanstd(classification_accuracy(:));
    classification_accuracy_stderr(t2) = classification_accuracy_std(t2)./sqrt(iteration);
end

    disp(['completing process...'])    
    cd(main_folder)    
    
% save data    
if flagTask == 1    
    foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Analysis/'];
    filename = ['spiral_classifier_performance.mat'];
    save([foldername,filename],'classification_accuracy','classification_accuracy_avg','classification_accuracy_stderr')
elseif flagTask == 2
    foldername = [main_folder,'/Sample Data/Language Task Additional 100 sub/Analysis/'];
    filename = ['spiral_classifier_performance.mat'];
    save([foldername,filename],'classification_accuracy','classification_accuracy_avg','classification_accuracy_stderr')
end

[classification_accuracy_avg_max t_max]= nanmax(classification_accuracy_avg(:));
classification_accuracy_stderr_max = classification_accuracy_stderr(t_max);
classification_accuracy_max = classification_accuracy(:,t_max);

figure()
bar(classification_accuracy_avg_max,'g')
% somenames={'original 100'; 'additional 100'; 'amplitude based'; } 
somenames={'original 100'} % language task
set(gca,'xticklabel',somenames)
hold on
er = errorbar(1:1,classification_accuracy_avg_max,classification_accuracy_stderr_max,classification_accuracy_stderr_max);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
hold on
    scatter(ones(size(classification_accuracy_max(:))).*(1+(rand(size(classification_accuracy_max(:)))-(3.5-3))/3),classification_accuracy_max(:),'r','filled')
%     scatter(ones(size(classification_accuracy_language(:))).*(1+(rand(size(classification_accuracy_language(:)))-(3.5-3))/3),classification_accuracy_language(:),'r','filled')
%     scatter(ones(size(classification_accuracy_language_add100(:))).*(1+(rand(size(classification_accuracy_language_add100(:)))-(3.5-3*2))/3),classification_accuracy_language_add100(:),'r','filled')
%     scatter(ones(size(classification_accuracy_language_raw_amp(:))).*(1+(rand(size(classification_accuracy_language_raw_amp(:)))-(3.5-3*3))/3),classification_accuracy_language_raw_amp(:),'r','filled')
hold off
% title(['working memory task, spiral classification performance'])
title(['spiral classification performance'])

end

