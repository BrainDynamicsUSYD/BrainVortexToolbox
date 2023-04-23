function [score,latent] = PCA(flagSur,hemisphere,main_folder,No_of_Subject,flagTask,listen_or_answer,motor_or_PCC);

%% PCA analysis of raw amplitude map, mathlisten 1-20

    % load task-evoked unfiltered trial averaged fMRI signal file
    foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Analysis'];
    cd(foldername)
    filename = ['task_evoked_unfiltered_fMRI_signal.mat'];
    load(filename) 
    
    % min-max normalization
    
    data_math_listen_avg_max = nanmax(data_math_listen_avg,[],3);
    data_math_listen_avg_min = nanmin(data_math_listen_avg,[],3);
    data_math_listen_avg_minmaxnorm = (data_math_listen_avg-data_math_listen_avg_min)./(data_math_listen_avg_max-data_math_listen_avg_min);
    
    data_math_answer_avg_max = nanmax(data_math_answer_avg,[],3);
    data_math_answer_avg_min = nanmin(data_math_answer_avg,[],3);
    data_math_answer_avg_minmaxnorm = (data_math_answer_avg-data_math_answer_avg_min)./(data_math_answer_avg_max-data_math_answer_avg_min);
    

%%
if hemisphere == 1
    if listen_or_answer == 1
       if motor_or_PCC == 1
        temp1 = permute(data_math_listen_avg_minmaxnorm(90:120,110:140,:),[3,1,2]); % left-motor/premotor
       elseif motor_or_PCC == 2
        temp1 = permute(data_math_listen_avg_minmaxnorm(128:158,175:205,:),[3,1,2]); % left-PCC-listen
       end
    elseif listen_or_answer == 2
    temp1 = permute(data_math_answer_avg_minmaxnorm(90:120,110:140,:),[3,1,2]); % left-motor/premotor
    end 
elseif hemisphere == 2
    if listen_or_answer == 1
        if motor_or_PCC == 1
           temp1 = permute(data_math_listen_avg_minmaxnorm(105:135,120:150,:),[3,1,2]); % right-listen Mor
        elseif motor_or_PCC == 2
           temp1 = permute(data_math_listen_avg_minmaxnorm(128:158,65:95,:),[3,1,2]); % right-listen_PCC
        end
    elseif listen_or_answer == 2
    temp1 = permute(data_math_answer_avg_minmaxnorm(95:125,120:150,:),[3,1,2]); % right-listen PCC
    end        
end
temp2 = reshape(temp1,20,[]);

[coeff,score,latent,tsquared,explained,mu]  = pca(temp2);

figure();
hold on
for iTime = 1:1 %size(score,1)
if hemisphere == 1
    if listen_or_answer == 1
        p = plot3(score(iTime:iTime+19,2),score(iTime:iTime+19,1),score(iTime:iTime+19,3),'b')
        p.LineWidth = 2;
        scatter3(score(1,2),score(1,1),score(1,3),35,'b','filled')
    elseif listen_or_answer == 2
        p = plot3(score(iTime:iTime+19,2),score(iTime:iTime+19,1),score(iTime:iTime+19,3),'b-')
        p.LineWidth = 2;
        scatter3(score(1,2),score(1,1),score(1,3),35,'b','filled')
    end
elseif hemisphere == 2
    if listen_or_answer == 1    
        p = plot3(score(iTime:iTime+19,2),score(iTime:iTime+19,1),score(iTime:iTime+19,3),'r')    
        scatter3(score(1,2),score(1,1),score(1,3),35,'r','filled')
        p.LineWidth = 2;
    elseif listen_or_answer == 2
        p = plot3(score(iTime:iTime+19,2),score(iTime:iTime+19,1),score(iTime:iTime+19,3),'r-')   
        p.LineWidth = 2;
        scatter3(score(1,2),score(1,1),score(1,3),35,'r','filled')
    end
end
    xlim([-15,15])
    ylim([-15,15])
% title(['task-evoked raw amplitude, listen 1-20, RIGHT PCC, PCA, PC1=',num2str(latent(1)),'% variance, PC2=',num2str(latent(2)),'%'])
title(['task-evoked raw amplitude, math listen 1-20, LEFT M1-PMd, PCA, PC1=',num2str(latent(1)),'% variance, PC2=',num2str(latent(2)),'%'])
end


end