function [task_type_name] = TaskLabel_Extract(flagRest,main_folder,No_of_Subject,flagTask)

%%

if flagTask == 1  % language task: original 100 sub
    TaskLabel_AllSubject_language_orig = [];
    for subject = 1:No_of_Subject
        foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Raw Data'];
        cd(foldername)
        name = dir([pwd]) ;
        foldername = [name(subject+2).folder,'/',name(subject+2).name,'/EVs'];
        cd(foldername)          
        math_listen = load('present_math.txt') ;
        math_question = load('question_math.txt') ;
        math_answer = load('response_math.txt') ;
        story_listen = load('present_story.txt') ;
        story_question = load('question_story.txt') ;
        story_answer = load('response_story.txt') ;  
        math_listen_restructure = [math_listen(:,1) math_listen(:,1)+math_listen(:,2) 4.*ones(size(math_listen,1),1)];
        math_question_restructure = [math_question(:,1) math_question(:,1)+math_question(:,2) 5.*ones(size(math_question,1),1)];
        math_answer_restructure = [math_answer(:,1) math_answer(:,1)+math_answer(:,2) 6.*ones(size(math_answer,1),1)];
        story_listen_restructure = [story_listen(:,1) story_listen(:,1)+story_listen(:,2) 1.*ones(size(story_listen,1),1)];
        story_question_restructure = [story_question(:,1) story_question(:,1)+story_question(:,2) 2.*ones(size(story_question,1),1)];
        story_answer_restructure = [story_answer(:,1) story_answer(:,1)+story_answer(:,2) 3.*ones(size(story_answer,1),1)];    
        
        TaskLabel_1sub = [math_listen_restructure;math_question_restructure;math_answer_restructure...
            ;story_listen_restructure;story_question_restructure;story_answer_restructure];
        TaskLabel_1sub_timestep = round([TaskLabel_1sub(:,1)./0.72 TaskLabel_1sub(:,2)./0.72 TaskLabel_1sub(:,3)]);
        [temp1,idx_sort] = sort(TaskLabel_1sub_timestep(:,1));
        TaskLabel_1sub_timestep_sort = TaskLabel_1sub_timestep(idx_sort,:);
        
        TaskLabel_AllSubject_language_orig{subject} = TaskLabel_1sub_timestep_sort;
    end
    task_type_name = {'present story','question story','answer story','present math','question math','answer math','cue','null'};
    foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Task Label/'];
    save([foldername,'LanguageOrigTaskLabelAllSubject.mat'],'TaskLabel_AllSubject_language_orig','task_type_name');

    
elseif flagTask == 2 % language task: additional 100 sub
    TaskLabel_AllSubject_language_add = [];
    for subject = 1:No_of_Subject

        foldername = [main_folder,'/Sample Data/Language Task Additional 100 sub/Raw Data'];
        cd(foldername)    
        name = dir([pwd]) ; 
        foldername = [name(subject+2).folder,'/',name(subject+2).name,'/EVs'];
        cd(foldername)    
        
        math_listen = load('present_math.txt') ;
        math_question = load('question_math.txt') ;
        math_answer = load('response_math.txt') ;
        story_listen = load('present_story.txt') ;
        story_question = load('question_story.txt') ;
        story_answer = load('response_story.txt') ;  
        math_listen_restructure = [math_listen(:,1) math_listen(:,1)+math_listen(:,2) 4.*ones(size(math_listen,1),1)];
        math_question_restructure = [math_question(:,1) math_question(:,1)+math_question(:,2) 5.*ones(size(math_question,1),1)];
        math_answer_restructure = [math_answer(:,1) math_answer(:,1)+math_answer(:,2) 6.*ones(size(math_answer,1),1)];
        story_listen_restructure = [story_listen(:,1) story_listen(:,1)+story_listen(:,2) 1.*ones(size(story_listen,1),1)];
        story_question_restructure = [story_question(:,1) story_question(:,1)+story_question(:,2) 2.*ones(size(story_question,1),1)];
        story_answer_restructure = [story_answer(:,1) story_answer(:,1)+story_answer(:,2) 3.*ones(size(story_answer,1),1)];    
        
        TaskLabel_1sub = [math_listen_restructure;math_question_restructure;math_answer_restructure...
            ;story_listen_restructure;story_question_restructure;story_answer_restructure];
        TaskLabel_1sub_timestep = round([TaskLabel_1sub(:,1)./0.72 TaskLabel_1sub(:,2)./0.72 TaskLabel_1sub(:,3)]);
        [temp1,idx_sort] = sort(TaskLabel_1sub_timestep(:,1));
        TaskLabel_1sub_timestep_sort = TaskLabel_1sub_timestep(idx_sort,:);
        
        TaskLabel_AllSubject_language_add{subject} = TaskLabel_1sub_timestep_sort;
    end
    task_type_name = {'present story','question story','answer story','present math','question math','answer math','cue','null'};
    foldername = [main_folder,'/Sample Data/Language Task Additional 100 sub/Task Label/'];
    save([foldername,'LanguageAddTaskLabelAllSubject.mat'],'TaskLabel_AllSubject_language_add','task_type_name');

    
elseif flagTask == 3 % working memory task
 
    TaskLabel_AllSubject_WM = [];
    for subject = 1:No_of_Subject
        foldername = [main_folder,'/Sample Data/Working Memory Task/Raw Data'];
        cd(foldername)
        name = dir([pwd]) ; 
        foldername = [name(subject+2).folder,'/',name(subject+2).name,'/EVs'];
        cd(foldername)
        
        bk0_body = load('0bk_body.txt') ;  
        bk0_body(:,4) = 1;
        bk0_cor = load('0bk_cor.txt');
        bk0_cor(:,4) = -1;
        bk0_err = load('0bk_err.txt');
        bk0_err(:,4) = -2;
        bk0_faces = load('0bk_faces.txt');
        bk0_faces(:,4) = 2;
        bk0_nlr = load('0bk_nlr.txt');
        bk0_nlr(:,4) = 9;
        bk0_places = load('0bk_places.txt');
        bk0_places(:,4) = 3;
        bk0_tools = load('0bk_tools.txt');
        bk0_tools(:,4) = 4;
        bk2_body = load('2bk_body.txt');
        bk2_body(:,4) = 5;
        bk2_cor = load('2bk_cor.txt');
        bk2_cor(:,4) = -3;
        bk2_err = load('2bk_err.txt');
        bk2_err(:,4) = -4;        
        bk2_faces = load('2bk_faces.txt');
        bk2_faces(:,4) = 6;  
        bk2_nlr = load('2bk_nlr.txt');
        bk2_nlr(:,4) = 10;
        bk2_places = load('2bk_places.txt');
        bk2_places(:,4) = 7;
        bk2_tools = load('2bk_tools.txt');
        bk2_tools(:,4) = 8;
        all_bk_cor = load('all_bk_cor.txt');
        all_bk_cor(:,4) = 11;
        all_bk_err = load('all_bk_err.txt');
        all_bk_err(:,4) = 12;

       TaskLabel = [bk0_body; bk0_cor; bk0_err; bk0_faces; bk0_nlr;bk0_places;bk0_tools;bk2_body;bk2_cor;bk2_err;bk2_faces;bk2_nlr;bk2_places;bk2_tools;all_bk_cor;all_bk_err];
       [a idx] = sort(TaskLabel(:,1));
       TaskLabel_sort = TaskLabel(idx,:);

       for label_i = 1:size(TaskLabel_sort,1)
           temp1 = TaskLabel_sort(label_i,4);
           if temp1 == 11 || temp1 == 12 || temp1 == 10
               continue
           end
           if temp1>0 
               temp2 = temp1;
           elseif temp1<0 
               temp2 = temp3; % apply the previous task label (task-block label) if is 0bk_cor or 0bk_err
           end
           if temp2~=9 && temp2~= 10
              temp3 = temp2; 
           end
           TaskLabel_sort(label_i,5) =temp3;


       end
       temp1 = TaskLabel_sort(:,5);
       idx = find(temp1>0); 
       TaskLabel_sort2 = TaskLabel_sort(idx,:); % remove 0 values/all-back cor/err tasks
       temp1 = TaskLabel_sort2(:,4);
       idx = find(temp1<0);
       TaskLabel_sort3 = TaskLabel_sort2(idx,:); % remove non-answer tasks-i.e., taskbloack onsets
       % restructure task label format
       onset = TaskLabel_sort3(:,1); % onset time
       offset = TaskLabel_sort3(:,1)+2.5; % offset time
       answer_cor_err = TaskLabel_sort3(:,4); 
       answer_cor_err(answer_cor_err==-1) = 1; % correct
       answer_cor_err(answer_cor_err==-2) = 0; % error
       answer_cor_err(answer_cor_err==-3) = 1; % correct
       answer_cor_err(answer_cor_err==-4) = 0; % error
       % task type:
       % 1=bk0_body;2=bk0_faces;3=bk0_places;4=bk0_tools;5=bk2_body;6=bk2_faces;7bk2_places=;8=bk2_tools
       task_type_name = ['1=bk0_body;2=bk0_faces;3=bk0_places;4=bk0_tools;5=bk2_body;6=bk2_faces;7bk2_places=;8=bk2_tools'];
       task_type = TaskLabel_sort3(:,5); 
       TaskLabel_restruture = [onset./0.72 offset./0.72 answer_cor_err task_type];
       TaskLabel_AllSubject_WM{subject} = TaskLabel_restruture;

    end
    % TaskLabel_AllSubject(TaskLabel_AllSubject==-1) = 1; % correct
    % TaskLabel_AllSubject(TaskLabel_AllSubject==-2) = 0; % error
    % TaskLabel_AllSubject(TaskLabel_AllSubject==-3) = 1; % correct
    % TaskLabel_AllSubject(TaskLabel_AllSubject==-4) = 0; % error

    foldername = [main_folder,'/Sample Data/Working Memory Task/Task Label/'];
    save([foldername,'WorkingMemoryTaskLabelAllSubject.mat'],'TaskLabel_AllSubject_WM','task_type_name');

end

cd(main_folder)

end