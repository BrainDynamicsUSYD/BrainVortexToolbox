function [contrast_density_parcelavg, contrast_density_parcel_stderr,contrast] = vortex_contrast_distribution_7networks(No_of_Subject,listen_or_answer,vortex_distribution_math,vortex_distribution_story)
disp(['initiating process...'])    
disp(['calculating vortex contrast (math vs story) distribution in 7 functional networks...'])    

%% contrast distribution in 7 functional networks

% calculate the subject-averaged vortex distribution: math tasks only 
vortex_distribution_story_matrix = [];
vortex_distribution_math_matrix = [];
empty_trial_vortex_math = [];
empty_trial_vortex_story = [];
vortex_distribution_math_matrix_subjectavg = [];
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
    vortex_distribution_math_matrix_subjectavg(:,:,:,subject) = nanmean(vortex_distribution_math_matrix,4);
end

% calculate the subject-averaged vortex distribution: story tasks only 
vortex_distribution_story_matrix_subjectavg = [];
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
    vortex_distribution_story_matrix_subjectavg(:,:,:,subject) = nanmean(vortex_distribution_math_matrix,4);

end

%% visualization: vortex contrast distribution in 7 functional networks
disp(['calculating the subject-averaged contrast density for each functional network...'])    

% load parcellation template of 7 functional networks
load('parcellation_template7.mat')

contrast = [];
for subject = 1:No_of_Subject
    contrast(:,:,:,subject) = abs(nanmean(vortex_distribution_story_matrix_subjectavg(:,:,:,subject),3) - nanmean(vortex_distribution_math_matrix_subjectavg(:,:,:,subject),3));
end

% calculate the mean contrast density for each of 7 functional networks (parcellation) 
parcellations = [1,2,3,4,5,6,7];
for p = 1:size(parcellations,2)
    for subject = 1:No_of_Subject 
        temp2_avg_norm_filt_1sub = contrast(:,:,:,subject);
        temp1_ind = temp2_avg_norm_filt_1sub(parcellation_template7 ==parcellations(p));   
        contrast_density_parcel(subject) = nansum(temp1_ind(:))./size(temp1_ind,1); % contrast density of a parcellation for each subject
    end
        contrast_density_parcelavg(p) = nanmean(contrast_density_parcel(:)); % mean contrast density across subjects 
        contrast_density_parcel_stderr(p) = nanstd(contrast_density_parcel(:))./sqrt(No_of_Subject); % standard error
end

% set error bar for standard error
errlow = - contrast_density_parcel_stderr;
errhigh =  contrast_density_parcel_stderr;
 
x = [1:7]; % x-axis
figure()
bar(x,contrast_density_parcelavg)
plotName = [{'VIS'},{'SMN'},{'AUD'},{'CON'},{'DAN'},{'FPN'},{'DMN'}] ;
set(gca,'xticklabel',plotName)
ylabel(['Contrast Density'])
hold on
er = errorbar(x,contrast_density_parcelavg,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
hold off
if listen_or_answer == 1
    title(['population-averaged vortex contrast density in 7 functional networks, Story vs math listening'])
elseif listen_or_answer == 2
    title(['population-averaged vortex contrast density in 7 functional networks, Story vs math answering'])
end
disp(['completing process...'])    

end