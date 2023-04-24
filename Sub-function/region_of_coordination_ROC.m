function [region_of_coordination] = region_of_coordination_ROC(flagSur,hemisphere,main_folder,No_of_Subject,flagTask);


%% task specific phase gradient field: math present 6-10 vs math answer 1-5: angle dif
% load task label         
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
 

cd(main_folder)   
load('parcellation_template7.mat')

Vx_math_listen_accu = [];
Vy_math_listen_accu = [];
trial_count_math_listen = 0;
trial_count_math_answer = 0;
Vx_math_answer_accu = [];
Vy_math_answer_accu = [];
for subject = 1:No_of_Subject
% load spatiotemporal bandpass filtered fMRI signal file
foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
cd(foldername)
filename = ['Preprocessed_spatiotemporalbandpass_data_language_task_orig100_LEFT_sub',num2str(subject),'.mat'];
load(filename) 
    
sigBPass = permute(DataOut,[2,3,4,1]);

% phase field
phaseSig = nan(size(sigBPass));
for irow = 1:size(sigBPass,1)
    for icol = 1:size(sigBPass,2)
        temp1 = sigBPass(irow,icol,:);
        phaseSig(irow,icol,:) = angle(hilbert(temp1(:)));
    end
end

% phase gradient (vector) field
Vx_flowmap_norm_phase = [];
Vy_flowmap_norm_phase = [];
vPhaseX = zeros(size(phaseSig)) ;
vPhaseY = zeros(size(phaseSig)) ;
cd(main_folder)
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

% task-specific phase vector field

    t_duration = 316;
    window_size = 20;
    temp1 = fullTime_allsubject{subject}; % task label
    temp1_time = temp1(:,3);
    % time segments for 4 tasks
    count = find(temp1_time==1);
    start_end_time_story_listen = temp1(count,:);
    count = find(temp1_time==3);
    start_end_time_story_answer= temp1(count,:);        
    count = find(temp1_time==4);
    start_end_time_math_listen = temp1(count,:);
    count = find(temp1_time==6);
    start_end_time_math_answer = temp1(count,:);        

    if nansum(Vx_flowmap_norm_phase(:))== 0
        continue
    end
    for trial = 1:size(start_end_time_math_listen,1);
        temp2_start = start_end_time_math_listen(trial,1)+ 1;
        temp2_end = temp2_start + window_size -1;
        if temp2_end > t_duration
            break
        end
        trial_count_math_listen = trial_count_math_listen + 1;            
        Vx_math_listen_accu(:,:,:,trial_count_math_listen) = Vx_flowmap_norm_phase(:,:,temp2_start:temp2_end);
        Vy_math_listen_accu(:,:,:,trial_count_math_listen) = Vy_flowmap_norm_phase(:,:,temp2_start:temp2_end);
    end  

    for trial = 1:size(start_end_time_math_answer,1);
        temp2_start = start_end_time_math_answer(trial,1)+ 1;
        temp2_end = temp2_start + window_size -1;
        if temp2_end > t_duration
            break
        end
        trial_count_math_answer = trial_count_math_answer + 1;              
        Vx_math_answer_accu(:,:,:,trial_count_math_answer) = Vx_flowmap_norm_phase(:,:,temp2_start:temp2_end);
        Vy_math_answer_accu(:,:,:,trial_count_math_answer) = Vy_flowmap_norm_phase(:,:,temp2_start:temp2_end);    end   
end

Vx_math_listen_accu_avg = nanmean(Vx_math_listen_accu,4);
Vy_math_listen_accu_avg = nanmean(Vy_math_listen_accu,4);
Vx_math_answer_accu_avg = nanmean(Vx_math_answer_accu,4);
Vy_math_answer_accu_avg = nanmean(Vy_math_answer_accu,4);

Vx_math_listen_accu_avg_6to10 = nanmean(Vx_math_listen_accu_avg(:,:,6:10),3);
Vy_math_listen_accu_avg_6to10 = nanmean(Vy_math_listen_accu_avg(:,:,6:10),3);
Vx_math_answer_accu_avg_1to5 = nanmean(Vx_math_answer_accu_avg(:,:,1:5),3);
Vy_math_answer_accu_avg_1to5 = nanmean(Vy_math_answer_accu_avg(:,:,1:5),3);

Vxy_math_listen_accu_avg_angle = angle(Vx_math_listen_accu_avg_6to10 + i.*Vy_math_listen_accu_avg_6to10);
Vxy_math_answer_accu_avg_angle = angle(Vx_math_answer_accu_avg_1to5 + i.*Vy_math_answer_accu_avg_1to5);

%% save data

save_folder = [main_folder,'/Sample Data/Language Task Original 100 sub/Analysis/'];
save([save_folder,'task_specific_phase_vector_field.mat'],'Vx_math_listen_accu','Vy_math_listen_accu','Vx_math_answer_accu','Vy_math_answer_accu')   ;   

%% region of coordination
angle_dif = anglesubtract(Vxy_math_listen_accu_avg_angle, Vxy_math_answer_accu_avg_angle) ;

angle_dif_abs = abs(angle_dif);
angle_dif_abs_filt = angle_dif_abs;
angle_dif_abs_filt(angle_dif_abs_filt<2/4*pi) = 0;
angle_dif_abs_filt(isnan(angle_dif_abs_filt)) = 0;
angle_dif_abs_filt(:,:,2) = 0;

params.minPattTime = 1;
params.minPattSize = 18; % 36mm
[WCentroids,absoluteTime,instantTotalPower,pattSize,patternIdx] = pattDetection_v4(angle_dif_abs_filt,angle_dif_abs_filt,params,0,'CenterOfMass_Amp');

angle_dif_abs_filt2 = nan(176,251);
for ipatt = 1:size(patternIdx,2)
    idx = patternIdx{ipatt};
    [y,x] = ind2sub([176,251,2],idx);
    for i2 = 1:size(y,1)
        angle_dif_abs_filt2(y(i2),x(i2)) = angle_dif_abs(y(i2),x(i2));
    end
% % 
% idx = patternIdx{1};
% [y,x] = ind2sub([176,251,2],idx);
% for i2 = 1:size(y,1)
%     angle_dif_abs_filt2(y(i2),x(i2)) = angle_dif_abs(y(i2),x(i2));
% end

end
%% angle difference map visualziation

region_of_coordination = angle_dif_abs_filt2;

xlimit_low = 130;
xlimit_high = 220;
ylimit_low = 20;
ylimit_high = 160;

figure(3)
subplot(1,2,1)
pcolor(angle_dif_abs_filt2)
shading interp
colormap autumn
colorbar
title(['math listen 6-10, trial avg phase vector and ROC'])
 xlim([xlimit_low, xlimit_high]);
caxis([0.5*pi,pi])
hold on
for parcellation_ID = 1:22
parcellation_template_1par = parcellation_template7;
parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
%     parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;
B = bwboundaries(parcellation_template_1par,'noholes');
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2),boundary(:,1),'-','linewidth',2,'color',[0,0,0])

end
end  
 set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
 hold off
 hold on
h = streamslice([1:251],[1:176],Vx_math_listen_accu_avg_6to10,Vy_math_listen_accu_avg_6to10,6,'k');
set( h, 'Color', [0 0 0] )
set(h,'linewidth',1.5)
hold off

subplot(1,2,2)
pcolor(angle_dif_abs_filt2)
shading interp
colormap autumn
colorbar
title(['math answer 1-5, trial avg phase vector and ROC'])
 xlim([xlimit_low, xlimit_high]);
caxis([0.5*pi,pi])
hold on
for parcellation_ID = 1:22
parcellation_template_1par = parcellation_template7;
parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
%     parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;
B = bwboundaries(parcellation_template_1par,'noholes');
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2),boundary(:,1),'-','linewidth',2,'color',[0,0,0])

end
end  
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
 hold off
 hold on
h = streamslice([1:251],[1:176],Vx_math_answer_accu_avg_1to5,Vy_math_answer_accu_avg_1to5,6,'k');
set( h, 'Color', [0 0 0] )
set(h,'linewidth',1.5)
hold off


end
