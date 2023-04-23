
function [spiral_distribution_math_listen_mid_avg,spiral_distribution_math_answer_end_avg,spiral_distribution_story_listen_mid_avg,spiral_distribution_story_answer_end_avg] = task_specific_spiral_distribution(flagRest,flagSur,hemisphere,main_folder,No_of_Subject,flagTask)

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

%% Task specific spiral centre distribution of each subject

if flagTask == 1 % language task: original 100 sub ONLY (for demonstration)
    for flagSur = 0:1 % both real and surrogate data are needed for further analysis
    % define parameters
    session_duration = 20;
    disp(['extracting task specific spiral centre distributions...'])    
    spiral_distribution_math_listen = [];
    spiral_distribution_math_answer = [];
    spiral_distribution_story_listen = [];
    spiral_distribution_story_answer = [];
    spiral_distribution_math_question = [];  
    spiral_distribution_story_question = [];  
    for subject = 1:No_of_Subject 
        foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Spiral Detected'];
        cd(foldername)
        filename = ['Spiral_detected_surfilt_language_task_orig100_LEFT_sub',num2str(subject),'.mat'];
        load(filename)
        if flagSur == 0
            spiral_filt_nega = spiral_filt_nega_real_95perc_extend;
            spiral_filt_pos = spiral_filt_pos_real_95perc_extend;
        elseif flagSur == 1
            spiral_filt_nega = spiral_filt_nega_sur;
            spiral_filt_pos = spiral_filt_pos_sur;        
        end

        if flagRest == 0
            temp1_accu_1subject_spiral = zeros(175,251,315);
        elseif flagRest == 1
            temp1_accu_1subject_spiral = zeros(175,251,1199);        
        end
        % clockwise spirals
        for ipatt = 1:size(spiral_filt_nega,1)
            for t = 1:size(spiral_filt_nega,2)
                temp1 = full(spiral_filt_nega{ipatt,t});
                if nansum(temp1(:))~=0
                   temp1_accu_1subject_spiral(:,:,t) = temp1_accu_1subject_spiral(:,:,t) + temp1;
                end
            end
        end
        % anticlockwise spirals
        for ipatt = 1:size(spiral_filt_pos,1)
            for t = 1:size(spiral_filt_pos,2)
                temp1 = full(spiral_filt_pos{ipatt,t});
                if nansum(temp1(:))~=0
                   temp1_accu_1subject_spiral(:,:,t) = temp1_accu_1subject_spiral(:,:,t) + temp1;
                end
            end
        end         

        t_duration = size(temp1_accu_1subject_spiral,3)-session_duration+1; % last time point to be used for analysis
        % extract task label for each subject
        temp1 = fullTime_allsubject{subject};
        temp1_time = temp1(:,3);
        % find time points for math tasks (listening or answering)
        count = find(temp1_time==4); % 4 for math listening tasks
        start_end_time_math_listen = temp1(count,:);
        count = find(temp1_time==5); % 4 for math listening tasks
        start_end_time_math_question = temp1(count,:);                
        count = find(temp1_time==6); % 6 for math answering tasks
        start_end_time_math_answer = temp1(count,:);

        % time points for story tasks (listening or answering)
        count = find(temp1_time==1); % label-1 for story listening tasks
        start_end_time_story_listen = temp1(count,:);
        count = find(temp1_time==2); % label-1 for story listening tasks
        start_end_time_story_question = temp1(count,:);                
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
                temp2_window = [temp2_start:temp2_start+session_duration-1];
                if temp2_start<1
                    continue
                end
                if temp2_end > t_duration
                    break
                end
                spiral_distribution_math_listen{trial,subject} = temp1_accu_1subject_spiral(:,:,t_select);
                if temp2_window(end) < t_duration
                    spiral_distribution_math_listen_1to20{trial,subject} = temp1_accu_1subject_spiral(:,:,temp2_window);
                end
            end
            for trial = 1:size(start_end_time_math_question,1);
                temp2_start = start_end_time_math_question(trial,1);
                temp2_end = start_end_time_math_question(trial,2);;
                temp2_mid = round((temp2_end-temp2_start)./2)+temp2_start;
                temp2_25perc = round((temp2_end-temp2_start).*0.25)+temp2_start;
                temp2_75perc = round((temp2_end-temp2_start).*0.75)+temp2_start;
                t_select = [temp2_start temp2_25perc temp2_mid temp2_75perc temp2_end];
                temp2_window = [temp2_start:temp2_start+session_duration-1];            
                if temp2_start<1
                    continue
                end
                if temp2_end > t_duration
                    break
                end
                spiral_distribution_math_question{trial,subject} = temp1_accu_1subject_spiral(:,:,t_select);
                if temp2_window(end) < t_duration
                    spiral_distribution_math_question_1to20{trial,subject} = temp1_accu_1subject_spiral(:,:,temp2_window);
                end
            end            
            for trial = 1:size(start_end_time_math_answer,1);
                temp2_start = start_end_time_math_answer(trial,1);                
                temp2_end = start_end_time_math_answer(trial,2);;
                temp2_mid = round((temp2_end-temp2_start)./2)+temp2_start;
                temp2_25perc = round((temp2_end-temp2_start).*0.25)+temp2_start;
                temp2_75perc = round((temp2_end-temp2_start).*0.75)+temp2_start;             
                t_select = [temp2_start temp2_25perc temp2_mid temp2_75perc temp2_end];
                temp2_window = [temp2_start:temp2_start+session_duration-1];            
                if temp2_start<1
                    continue
                end
                if temp2_end > t_duration 
                    break
                end
                spiral_distribution_math_answer{trial,subject} = temp1_accu_1subject_spiral(:,:,t_select);
                if temp2_window(end) < t_duration
                    spiral_distribution_math_answer_1to20{trial,subject} = temp1_accu_1subject_spiral(:,:,temp2_window);
                end
            end            
            % extract task specifc spiral distribution data: story tasks
            for trial = 1:size(start_end_time_story_listen,1);
                temp2_start = start_end_time_story_listen(trial,1);                
                temp2_end = start_end_time_story_listen(trial,2);;
                temp2_mid = round((temp2_end-temp2_start)./2)+temp2_start;
                temp2_25perc = round((temp2_end-temp2_start).*0.25)+temp2_start;
                temp2_75perc = round((temp2_end-temp2_start).*0.75)+temp2_start;
                t_select = [temp2_start temp2_25perc temp2_mid temp2_75perc temp2_end];
                temp2_window = [temp2_start:temp2_start+session_duration-1];            
                if temp2_start<1
                    continue
                end                
                if temp2_end > t_duration 
                    break
                end
                spiral_distribution_story_listen{trial,subject} = temp1_accu_1subject_spiral(:,:,t_select);
                if temp2_window(end) < t_duration            
                    spiral_distribution_story_listen_1to20{trial,subject} = temp1_accu_1subject_spiral(:,:,temp2_window);
                end
            end
            for trial = 1:size(start_end_time_story_question,1);
                temp2_start = start_end_time_story_question(trial,1);                
                temp2_end = start_end_time_story_question(trial,2);;
                temp2_mid = round((temp2_end-temp2_start)./2)+temp2_start;
                temp2_25perc = round((temp2_end-temp2_start).*0.25)+temp2_start;
                temp2_75perc = round((temp2_end-temp2_start).*0.75)+temp2_start;
                t_select = [temp2_start temp2_25perc temp2_mid temp2_75perc temp2_end];
                temp2_window = [temp2_start:temp2_start+session_duration-1];             
                if temp2_start<1
                    continue
                end                
                if temp2_end > t_duration 
                    break
                end
                spiral_distribution_story_question{trial,subject} = temp1_accu_1subject_spiral(:,:,t_select);
                if temp2_window(end) < t_duration               
                    spiral_distribution_story_question_1to20{trial,subject} = temp1_accu_1subject_spiral(:,:,temp2_window);
                end
            end            
            for trial = 1:size(start_end_time_story_answer,1);
                temp2_start = start_end_time_story_answer(trial,1);                
                temp2_end = start_end_time_story_answer(trial,2);;
                temp2_mid = round((temp2_end-temp2_start)./2)+temp2_start;
                temp2_25perc = round((temp2_end-temp2_start).*0.25)+temp2_start;
                temp2_75perc = round((temp2_end-temp2_start).*0.75)+temp2_start;
                t_select = [temp2_start temp2_25perc temp2_mid temp2_75perc temp2_end];
                temp2_window = [temp2_start:temp2_start+session_duration-1];             
                if temp2_start<1
                    continue
                end                
                if temp2_end > t_duration 
                    break
                end
                spiral_distribution_story_answer{trial,subject} = temp1_accu_1subject_spiral(:,:,t_select);
                if temp2_window(end) < t_duration                           
                    spiral_distribution_story_answer_1to20{trial,subject} = temp1_accu_1subject_spiral(:,:,temp2_window);
                end
            end            

    end


    no_of_trial = 0;
    for subject = 1:size(spiral_distribution_math_listen,2)
        for trial = 1:size(spiral_distribution_math_listen,1)
            temp1 = full(spiral_distribution_math_listen{trial,subject});
            if nansum(temp1(:))~=0
               no_of_trial = no_of_trial+ 1;
               spiral_distribution_math_listen_mid(:,:,no_of_trial) = temp1(:,:,3);           
            end       
        end
        spiral_distribution_math_listen_mid_avg = nanmean(spiral_distribution_math_listen_mid,3);
    end

    no_of_trial = 0;
    for subject = 1:size(spiral_distribution_math_answer,2)
        for trial = 1:size(spiral_distribution_math_answer,1)
            temp1 = full(spiral_distribution_math_answer{trial,subject});
            if nansum(temp1(:))~=0
               no_of_trial = no_of_trial+ 1;
               spiral_distribution_math_answer_end(:,:,no_of_trial) = temp1(:,:,5);           
            end       
        end
        spiral_distribution_math_answer_end_avg = nanmean(spiral_distribution_math_answer_end,3);
    end

    no_of_trial = 0;
    for subject = 1:size(spiral_distribution_story_listen,2)
        for trial = 1:size(spiral_distribution_story_listen,1)
            temp1 = full(spiral_distribution_story_listen{trial,subject});
            if nansum(temp1(:))~=0
               no_of_trial = no_of_trial+ 1;
               spiral_distribution_story_listen_mid(:,:,no_of_trial) = temp1(:,:,3);           
            end       
        end
        spiral_distribution_story_listen_mid_avg = nanmean(spiral_distribution_story_listen_mid,3);
    end

    no_of_trial = 0;
    for subject = 1:size(spiral_distribution_story_answer,2)
        for trial = 1:size(spiral_distribution_story_answer,1)
            temp1 = full(spiral_distribution_story_answer{trial,subject});
            if nansum(temp1(:))~=0
               no_of_trial = no_of_trial+ 1;
               spiral_distribution_story_answer_end(:,:,no_of_trial) = temp1(:,:,5);           
            end       
        end
        spiral_distribution_story_answer_end_avg = nanmean(spiral_distribution_story_answer_end,3);
    end
    if flagSur == 0
        foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Analysis/'];
        filename = ['task_specific_spiral_distribution.mat'];
        save([foldername,filename],'spiral_distribution_math_listen_mid','spiral_distribution_math_listen_mid_avg'...
            ,'spiral_distribution_math_answer_end_avg','spiral_distribution_math_answer_end'...
            ,'spiral_distribution_story_listen_mid','spiral_distribution_story_listen_mid_avg'...
            ,'spiral_distribution_story_answer_end','spiral_distribution_story_answer_end_avg'...
            ,'spiral_distribution_math_listen','spiral_distribution_math_answer','spiral_distribution_story_listen','spiral_distribution_story_answer'...
            ,'spiral_distribution_math_listen_1to20','spiral_distribution_math_answer_1to20','spiral_distribution_story_listen_1to20','spiral_distribution_story_answer_1to20');
    elseif flagSur == 1
        foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Analysis/'];
        filename = ['task_specific_spiral_distribution_sur.mat'];
        save([foldername,filename],'spiral_distribution_math_listen_mid','spiral_distribution_math_listen_mid_avg'...
            ,'spiral_distribution_math_answer_end_avg','spiral_distribution_math_answer_end'...
            ,'spiral_distribution_story_listen_mid','spiral_distribution_story_listen_mid_avg'...
            ,'spiral_distribution_story_answer_end','spiral_distribution_story_answer_end_avg'...
            ,'spiral_distribution_math_listen','spiral_distribution_math_answer','spiral_distribution_story_listen','spiral_distribution_story_answer'...
            ,'spiral_distribution_math_listen_1to20','spiral_distribution_math_answer_1to20','spiral_distribution_story_listen_1to20','spiral_distribution_story_answer_1to20');
    end
    cd(main_folder)


    if flagSur == 0
    figure()
    subplot(2,2,1)
    imagesc(spiral_distribution_math_listen_mid_avg)
    set(gca,'ydir','normal')
    colormap jet
    colorbar
    title(['math listen end'])
    subplot(2,2,2)
    imagesc(spiral_distribution_math_answer_end_avg)
    set(gca,'ydir','normal')
    colormap jet
    colorbar
    title(['math answer mid'])
    subplot(2,2,3)
    imagesc(spiral_distribution_story_listen_mid_avg)
    set(gca,'ydir','normal')
    colormap jet
    colorbar
    title(['story listen end'])
    subplot(2,2,4)
    imagesc(spiral_distribution_story_answer_end_avg)
    set(gca,'ydir','normal')
    colormap jet
    colorbar
    title(['story answer mid'])
    end

    end

end

end

