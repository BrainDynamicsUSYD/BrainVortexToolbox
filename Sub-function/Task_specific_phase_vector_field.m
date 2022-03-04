function [vx_math_listen_alltrials,vy_math_listen_alltrials,vx_math_answer_alltrials,vy_math_answer_alltrials] = Task_specific_phase_vector_field(No_of_Subject,listen_or_answer)
disp(['initiating...'])    
disp(['calculating task-specific phase vector field (math listening and answering)...'])    

%% task specific phase vector field (math listening vs answering tasks)

disp(['calculating phase vector field via phase map...'])    

vortex_filt_nega = [];
vortex_filt_pos = [];

% load preprocessed fMRI data files
for subject = 1:No_of_Subject
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

    Vx_flowmap_norm_phase(isnan(Vx_flowmap_norm_phase)) = 0;
    Vy_flowmap_norm_phase(isnan(Vy_flowmap_norm_phase)) = 0;

    clearvars vPhaseX  vPhaseY phaseSig 

    for time = 1:size(data_allsubject,3)
                temp1_Vx_sparse{time} = sparse(Vx_flowmap_norm_phase(:,:,time));
                temp1_Vy_sparse{time} = sparse(Vy_flowmap_norm_phase(:,:,time));
                temp1_Vy_sparse_accu{time,subject} = temp1_Vy_sparse{time};
                temp1_Vx_sparse_accu{time,subject} = temp1_Vx_sparse{time};
    end
    clearvars vPhaseX  vPhaseY phaseSig temp1_phase  temp1_Vx_sparse temp1_Vy_sparse  Vx_flowmap_norm_phase Vy_flowmap_norm_phase
end

%% load task label
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

%% calculate average phase vector field for each trial from each subject

disp(['extracting phase vector field of each task trial session...'])    

temp1_story_present_train_template_accu = [];
temp1_math_present_train_template_accu = [];
vx_story_avg = [];
vy_story_avg = [];
vx_math_avg = [];
vy_math_avg = [];

session_duration = 5; % duration of each task session used for analysis
t_duration = size(data_allsubject,3)-session_duration+1; % duration of each task session used for analysis
t1 = 6; % starting time of a trial session to be analyzed, 6 for listeing tasks
t2 = 1;% starting time of a trial session to be analyzed, 1 for answering tasks

    vx_math_listen = [];
    vy_math_listen = [];
    vx_math_answer = [];
    vy_math_answer = [];

    
for subject = 1:No_of_Subject
        temp1_Vx_sparse_accu_1subject = [];
        temp1_Vy_sparse_accu_1subject = [];
        % extract u,v coordinates across time for each subject
        for time = 1:t_duration
            temp1_Vx_sparse_accu_1subject(:,:,time) = full(temp1_Vx_sparse_accu{time,subject});
            temp1_Vy_sparse_accu_1subject(:,:,time) = full(temp1_Vy_sparse_accu{time,subject});
        end

        % extract task label and session time points for each subject        
        temp1 = fullTime_allsubject{subject};
        temp1_time = temp1(:,3);
        % find time points for math listening tasks 
        count = find(temp1_time==4);
        start_end_time_math_listen = temp1(count,:);

        % find time points for math answering tasks 
        count = find(temp1_time==6);
        start_end_time_math_answer = temp1(count,:);

        % extract task specifc phase vector field: math listening tasks
        for trial = 1:size(start_end_time_math_listen,1);
            temp2_start = start_end_time_math_listen(trial,1)+ t1;
            temp2_end = temp2_start+session_duration-1;
            if temp2_end > t_duration
                break
            end
            vx_math_listen{trial,subject} = temp1_Vx_sparse_accu_1subject(:,:,temp2_start:temp2_end);
            vy_math_listen{trial,subject} = temp1_Vy_sparse_accu_1subject(:,:,temp2_start:temp2_end);

        end

        % extract task specifc phase vector field: math answering tasks
        for trial = 1:size(start_end_time_math_answer,1);
            temp2_start = start_end_time_math_answer(trial,1)+t2;
            temp2_end = temp2_start + session_duration-1;
            if temp2_end > t_duration
                break
            end
            vx_math_answer{trial,subject} = temp1_Vx_sparse_accu_1subject(:,:,temp2_start:temp2_end);
            vy_math_answer{trial,subject} = temp1_Vy_sparse_accu_1subject(:,:,temp2_start:temp2_end);    
        end
end

    clearvars temp1_Vy_sparse_accu_1subject temp1_Vy_sparse_accu_1subject temp1_Vx_sparse_accu_1subject temp1_Vy_sparse_accu_1subject



%% rearrange task-specifc phase vector field data by trials  
disp(['rearranging phase vector field distribution by trials...'])    

    vx_math_listen_alltrials = [];
    vy_math_listen_alltrials = [];
    % math listening tasks only
    trialNo = 0;
    for trial = 1:size(vx_math_listen,1)
        for subject = 1:No_of_Subject
            temp1_vx_math_listen = vx_math_listen{trial,subject}; 
            temp1_vy_math_listen = vy_math_listen{trial,subject}; 
            if nansum(temp1_vx_math_listen(:)) ~= 0   
                trialNo = trialNo + 1;  
                vx_math_listen_alltrials(:,:,:,trialNo) = temp1_vx_math_listen;
                vy_math_listen_alltrials(:,:,:,trialNo) = temp1_vy_math_listen;
            end
        end
    end
    vy_math_listen_alltrials(vy_math_listen_alltrials==0) = nan;
    vx_math_listen_alltrials(vx_math_listen_alltrials==0) = nan;
    
    % math answering tasks only    
    trialNo = 0;
    for trial = 1:size(vx_math_answer,1)
        for subject = 1:No_of_Subject
            temp1_vx_math_answer = vx_math_answer{trial,subject}; 
            temp1_vy_math_answer = vy_math_answer{trial,subject}; 
            if nansum(temp1_vx_math_answer(:)) ~= 0 
                trialNo = trialNo + 1;  
                vx_math_answer_alltrials(:,:,:,trialNo) = temp1_vx_math_answer;
                vy_math_answer_alltrials(:,:,:,trialNo) = temp1_vy_math_answer;        
            end
        end
    end
    vy_math_answer_alltrials(vy_math_answer_alltrials==0) = nan;
    vx_math_answer_alltrials(vx_math_answer_alltrials==0) = nan;

    
%% calculate the trial averaged angle difference map between math listening and asnwering tasks

disp(['calculating trial-averaged phase vector field angle differences between tasks...'])    

% calculate the mean phase vector field for math listening and answering
% tasks
vx_math_listen_alltrials_avg = nanmean(nanmean(vx_math_listen_alltrials,4),3);
vy_math_listen_alltrials_avg = nanmean(nanmean(vy_math_listen_alltrials,4),3);

vx_math_answer_alltrials_avg = nanmean(nanmean(vx_math_answer_alltrials,4),3);
vy_math_answer_alltrials_avg = nanmean(nanmean(vy_math_answer_alltrials,4),3);

% standardize the vector lengths of all phase vectors as 1
vx_math_answer_alltrials_avg_01 = vx_math_answer_alltrials_avg./sqrt(vx_math_answer_alltrials_avg.*vx_math_answer_alltrials_avg + vy_math_answer_alltrials_avg.*vy_math_answer_alltrials_avg);
vy_math_answer_alltrials_avg_01 = vy_math_answer_alltrials_avg./sqrt(vx_math_answer_alltrials_avg.*vx_math_answer_alltrials_avg + vy_math_answer_alltrials_avg.*vy_math_answer_alltrials_avg);
vx_math_listen_alltrials_avg_01 = vx_math_listen_alltrials_avg./sqrt(vx_math_listen_alltrials_avg.*vx_math_listen_alltrials_avg + vy_math_listen_alltrials_avg.*vy_math_listen_alltrials_avg);
vy_math_listen_alltrials_avg_01 = vy_math_listen_alltrials_avg./sqrt(vx_math_listen_alltrials_avg.*vx_math_listen_alltrials_avg + vy_math_listen_alltrials_avg.*vy_math_listen_alltrials_avg);

% calculate the angles of the phase vectors for both math listening and
% answering tasks
angle_vxy_math_answer_avg = angle(vx_math_answer_alltrials_avg_01 + i.*vy_math_answer_alltrials_avg_01);
angle_vxy_math_listen_avg = angle(vx_math_listen_alltrials_avg_01 + i.*vy_math_listen_alltrials_avg_01);

% calculate angle difference map between phase vector fields of math
% listening and answering tasks
angle_dif = anglesubtract(angle_vxy_math_answer_avg, angle_vxy_math_listen_avg) ;

% isolate cortical regions with angle differences > 0.5pi and pattern size 
% > 50 voxels
angle_dif_abs = abs(angle_dif);
angle_dif_abs_filt = angle_dif_abs;
angle_dif_abs_filt(angle_dif_abs_filt<0.5*pi) = 0;
angle_dif_abs_filt(isnan(angle_dif_abs_filt)) = 0;

angle_dif_abs_filt(:,:,2) = 0; % set the 3rd dimension to 0 to enable pattDetection_v4.m to work
params.minPattTime = 1;
params.minPattSize = 50;  
[WCentroids,absoluteTime,instantTotalPower,pattSize,patternIdx] = pattDetection_v4(angle_dif_abs_filt,angle_dif_abs_filt,params,0,'CenterOfMass_Amp');

idx = patternIdx{2};
[y,x] = ind2sub([176,251,2],idx);
angle_dif_abs_filt2 = nan(176,251);
for i2 = 1:size(y,1)
    angle_dif_abs_filt2(y(i2),x(i2)) = angle_dif_abs(y(i2),x(i2));
end

idx = patternIdx{1};
[y,x] = ind2sub([176,251,2],idx);
for i2 = 1:size(y,1)
    angle_dif_abs_filt2(y(i2),x(i2)) = angle_dif_abs(y(i2),x(i2));
end         

%% visualization: plot the cortical region of significant (>0.5pi) angle differencse in a 2D heatmap 
disp(['visualzing trial-averaged phase vector field angle difference map...'])    

xlimit_low = 130;
xlimit_high = 220;

% 2D heat map of angle difference between math listening and answering
% tasks, overlayed with the borders of 7 parcellation template and the 
% streamlines from trial average phase vector field of 
% math listening tasks
figure()
pcolor(angle_dif_abs_filt2)
shading interp
colormap autumn
colorbar
title(['trial-averaged phase vector angle difference (math listening vs. answering, >0.5pi) + math listening streamlines'])
 xlim([xlimit_low, xlimit_high]);
caxis([0.5*pi,pi])
% the borders of 7 parcellation template
load('parcellation_template7.mat')
hold on
 for parcellation_ID = 1:22
    parcellation_template_1par = parcellation_template7;
    parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
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
% the streamlines from trial average phase vector field of 
% math listening tasks
h = streamslice([1:251],[1:176],vx_math_listen_alltrials_avg_01,vy_math_listen_alltrials_avg_01,6,'k');
set( h, 'Color', [0 0 0] )
set(h,'linewidth',1.5)
hold off

% 2D heat map of angle difference between math listening and answering
% tasks, overlayed with the borders of 7 parcellation template and the 
% streamlines from trial average phase vector field of 
% math answering tasks
figure()
pcolor(angle_dif_abs_filt2)
shading interp
colormap autumn
colorbar
title(['trial-averaged phase vector angle difference (math listening vs. answering, >0.5pi) + math answering streamlines'])
 xlim([xlimit_low, xlimit_high]);
caxis([0.5*pi,pi])
% the borders of 7 parcellation template
load('parcellation_template7.mat')
hold on
 for parcellation_ID = 1:22
    parcellation_template_1par = parcellation_template7;
    parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
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
% the streamlines from trial average phase vector field of 
% math listening tasks
h = streamslice([1:251],[1:176],vx_math_answer_alltrials_avg_01,vy_math_answer_alltrials_avg_01,6,'k');
set( h, 'Color', [0 0 0] )
set(h,'linewidth',1.5)
hold off

%% single sample voxel task specific phase vector field distribution (polarhistogram)
disp(['visualzing task-specific phase vector angle distribution of a sample voxel...'])    

% select the x,y coordinate of a sample voxel 
x = 170;
y = 56;

% calculate the angle of all phase vectors at that sample voxel in both
% math listeing and answering tasks
temp1_vxy_math_answer_angle = angle(vx_math_answer_alltrials + i.*vy_math_answer_alltrials);
temp1_vxy_math_listen_angle = angle(vx_math_listen_alltrials + i.*vy_math_listen_alltrials);

temp1_vxy_math_answer_angle_1voxel = temp1_vxy_math_answer_angle(y,x,:,:);
temp1_vxy_math_listen_angle_1voxel = temp1_vxy_math_listen_angle(y,x,:,:);
temp1_vxy_math_answer_angle_1voxel(temp1_vxy_math_answer_angle_1voxel==0) = nan;
temp1_vxy_math_listen_angle_1voxel(temp1_vxy_math_listen_angle_1voxel==0) = nan;
figure()
subplot(1,2,1)
polarhistogram(temp1_vxy_math_answer_angle_1voxel(:),12, 'Normalization','probability')
title(['single voxel (x/y=170/56) sample vector angle, math answering'])
subplot(1,2,2)
polarhistogram(temp1_vxy_math_listen_angle_1voxel(:),12, 'Normalization','probability')
title(['single voxel (x/y=170/56) sample vector angle, math listening'])

end