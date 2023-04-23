
function [temp1_phase,Vx_flowmap_norm_phase,Vy_flowmap_norm_phase] = task_evoked_unfiltered_fMRI_signal(flagSur,hemisphere,main_folder,No_of_Subject,flagTask,listen_or_answer);

%% load task label (language task)        

clearvars i
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
 

%%  Task trial averaged ampltiude with raw and bandpasssed data


    
min_duration = 1; % minimum duration of a spiral
data_math_listen_accu = [];
data_math_answer_accu = [];
trial_count_math_listen = 0;
trial_count_math_answer = 0;
spiral_template_posi_nega_accu = [];

for subject = 1:No_of_Subject 
    foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
    cd(foldername)
    if hemisphere == 1
        filename = ['Preprocessed_raw_data_language_task_orig100_LEFT_sub',num2str(subject),'.mat'];
    elseif hemisphere == 2
        filename = ['Preprocessed_raw_data_language_task_orig100_RIGHT_sub',num2str(subject),'.mat'];
    end
        load(filename) 
    signal_raw = DataOut;
    signal_raw_avg = nanmean(signal_raw,3);
    signal_raw_norm = signal_raw - signal_raw_avg;

    t_duration = 316;
    window_size = 20;
    temp1 = fullTime_allsubject{subject}; % task label
    temp1_time = temp1(:,3);
    data_math_listen = []; 
    data_math_answer = [];
    % time segments for 4 tasks
    count = find(temp1_time==1);
    start_end_time_story_listen = temp1(count,:);
    count = find(temp1_time==3);
    start_end_time_story_answer= temp1(count,:);        
    count = find(temp1_time==4);
    start_end_time_math_listen = temp1(count,:);
    count = find(temp1_time==6);
    start_end_time_math_answer = temp1(count,:);        

    signal_raw_norm_1sub = signal_raw_norm;
    if nansum(signal_raw_norm_1sub(:))== 0
        continue
    end
    for trial = 1:size(start_end_time_math_listen,1);
        temp2_start = start_end_time_math_listen(trial,1)+ 1;
        temp2_end = temp2_start + window_size -1;
        if temp2_end > t_duration
            break
        end
        trial_count_math_listen = trial_count_math_listen + 1;            
        data_math_listen_accu(:,:,:,trial_count_math_listen) = signal_raw_norm_1sub(:,:,temp2_start:temp2_end);
    end  

    for trial = 1:size(start_end_time_math_answer,1);
        temp2_start = start_end_time_math_answer(trial,1)+ 1;
        temp2_end = temp2_start + window_size -1;
        if temp2_end > t_duration
            break
        end
        trial_count_math_answer = trial_count_math_answer + 1;              
        data_math_answer_accu(:,:,:,trial_count_math_answer) = signal_raw_norm_1sub(:,:,temp2_start:temp2_end);
    end   
end
 
data_math_listen_avg = nanmean(data_math_listen_accu,4);
data_math_answer_avg = nanmean(data_math_answer_accu,4);

% trial averaged phase gradien field
temp1_phase = nan(175,251,window_size);

if listen_or_answer == 1
    % math listen = 1
    for irow = 1:175
        for icol = 1:251
            temp1 = data_math_answer_avg(irow,icol,:);
            temp1_phase(irow,icol,:) = angle(hilbert(temp1(:)));
        end
    end
    elseif listen_or_answer == 2
    % math answer = 2
    for irow = 1:175
        for icol = 1:251
            temp1 = data_math_answer_avg(irow,icol,:);
            temp1_phase(irow,icol,:) = angle(hilbert(temp1(:)));
        end
    end
end
% phase gradient velocity field: Rory's code: anglesubtract.m
Vx_flowmap_norm_phase = [];
Vy_flowmap_norm_phase = [];
cd(main_folder)

phaseSig = temp1_phase(:,:,:);
vPhaseX = zeros(size(phaseSig)) ;
vPhaseY = zeros(size(phaseSig)) ;
for iTime = 1:size(phaseSig,3)
     for iX = 1:size(phaseSig,1)
     vPhaseX(iX,2:end-1,iTime) = (anglesubtract(  phaseSig(iX,3:end,iTime),  phaseSig(iX,1:end-2,iTime)  ))/2 ;
     end
     for iY = 1:size(phaseSig,2)
     vPhaseY(2:end-1,iY,iTime) = (anglesubtract(phaseSig(3:end,iY,iTime), phaseSig(1:end-2,iY,iTime))  )/2 ;
     end
end
Vx_flowmap_norm_phase(:,:,:) = -vPhaseX./sqrt(vPhaseX.^2+vPhaseY.^2) ;
Vy_flowmap_norm_phase(:,:,:) = -vPhaseY./sqrt(vPhaseX.^2+vPhaseY.^2) ;

% Detect candidate brain spirals via curl value calcualted from phase vector field: anticlockwise spiral only  
% clockwise and anticlockwise spirals are detected treated seperately
disp(['detecting full-sized anticlockwise brain spirals via phase vector field...'])    

spiral_filt_pos = [];
spiral_filt_nega = [];

curlz = [];
for time = 1:size(temp1_phase,3)     
    temp1_vx = Vx_flowmap_norm_phase(:,:,time);
    temp1_vy = Vy_flowmap_norm_phase(:,:,time);
    [x,y] = meshgrid(1:size(temp1_vx,2),1:size(temp1_vx,1));
    [curlz(:,:,time) cav(:,:,time)]  = curl(x,y,temp1_vx,temp1_vy);  % curl value of phase vector field    
end      
cav_filt_pos = curlz;
cav_filt_pos(cav_filt_pos<1) = 0;      % detect only anticlockwise vorties
cav_filt_pos(isnan(cav_filt_pos)) = 0; % filter out data points with curl value < 1

% spiral core detection based on filtered curl values (requires function
% pattDetection_V4.m)

params.minPattTime = min_duration; % minimum duration of spiral (10 time steps)
params.minPattSize = 3;  % minimum size of spiral (3x3 matrix)

[WCentroids,absoluteTime,instantTotalPower,pattSize,patternIdx] = pattDetection_v4(cav_filt_pos,cav_filt_pos,params,0,'CenterOfMass_Amp');
clearvars instantTotalPower patternIdx pattSize temp1_vx temp1_vy   
 
% find the angle difference between center-originating vector (center of mass point to actual vector location) and phase vector

if hemisphere == 1 % left brain
load parcellation_template % load parcellation template to confirm the edge of flattened cortical map 
parcellation_template_01 = parcellation_template(1:175,1:251)./parcellation_template(1:175,1:251);
elseif hemisphere == 2 % righ brain
load parcellation_template22_RightBrain_subject1-100 % load right brain parcellation template 
parcellation_template_01 = parcellation_template22_RightBrain_100sub(1:175,:,1)./parcellation_template22_RightBrain_100sub(1:175,:,1);
end

for ipatt = 1:size(absoluteTime,2) % select 1 spiral pattern detected   
      for time_withinpatt = 1:size(absoluteTime{ipatt},2);  % select 1 time point within the duration of the detect spiral pattern
        clearvars i
        y_tempalte = [];
        x_tempalte = [];
        centerofmass_1patt_vxy_angle_dif = [];
        centerofmass_1patt_vxy_angle_dif_fil = [];

        temp1_time = absoluteTime{ipatt};
        time = temp1_time(time_withinpatt); 
        
        % find the angles of all the phase vectors within the defined 
        % vector at the defined time 
        temp1_vx = Vx_flowmap_norm_phase(:,:,time); 
        temp1_vy = Vy_flowmap_norm_phase(:,:,time);
        temp1_vxy_angle = angle(temp1_vx + 1i.*temp1_vy); 
        
        % centre of mass position of the defined vector at the definded
        % time
        temp1_centerofmass_1patt = WCentroids{ipatt};
        centerofmass_1patt = temp1_centerofmass_1patt(time_withinpatt,:);

        % calculate the angles of spiral-centre-originated vectors toward all 
        % voxel positions across the flattened cortex
        for irow = 1:size(temp1_phase,1) 
            for icol = 1:size(temp1_phase,2) 
                y_tempalte(irow,:) = irow;
                x_tempalte(:,icol) = icol;   
            end
        end

        centerofmass_1patt_vx = round(x_tempalte - centerofmass_1patt(1));
        centerofmass_1patt_vy = round(y_tempalte - centerofmass_1patt(2));
        centerofmass_1patt_vxy_angle = angle(centerofmass_1patt_vx + 1i.*centerofmass_1patt_vy);
        
        % lengths of spiral-centre-originated vectors/ distances
        centerofmass_1patt_vxy_abs = sqrt(centerofmass_1patt_vx.*centerofmass_1patt_vx+centerofmass_1patt_vy.*centerofmass_1patt_vy);
        centerofmass_1patt_vxy_abs(centerofmass_1patt_vxy_abs==0) = 1;
        
        % angle differences between ortex-centre-originated vectors and
        % phase vectors at the same voxel
        centerofmass_1patt_vxy_angle_dif = centerofmass_1patt_vxy_angle - temp1_vxy_angle;

        for irow = 1:size(Vx_flowmap_norm_phase,1)
            for icol = 1:size(Vx_flowmap_norm_phase,2)
                temp1 = centerofmass_1patt_vxy_angle_dif(irow,icol);
                if temp1> pi        % limit the angle differences between -pi and pi
                    centerofmass_1patt_vxy_angle_dif(irow,icol) = centerofmass_1patt_vxy_angle_dif(irow,icol) - 2*pi;
                elseif temp1< -pi   % limit the angle differences between -pi and pi
                    centerofmass_1patt_vxy_angle_dif(irow,icol) = centerofmass_1patt_vxy_angle_dif(irow,icol) + 2*pi;
                end
            end
        end
        centerofmass_1patt_vxy_angle_dif = abs(centerofmass_1patt_vxy_angle_dif); 

        % angle difference threshold: phase vectors with angle differences 
        % >135 or <45 degrees from the centre-orginated-vector are excluded 
        % from the spiral, ideally angle difference should be 90 degrees
        % for a perfect spiral
        
        angle_threshold_low = pi./2 - 0.7854; %***********%0.7854; %45degrees- %0.5236; % 30degrees below %0.2618; % 15 degrees below
        angle_threshold_high = pi./2 + 0.7854; %*********%0.7854; %45degrees+ %0.2536; % 30degrees above %0.2618; % 15 degrees above
        centerofmass_1patt_vxy_angle_dif_fil = centerofmass_1patt_vxy_angle_dif;
        centerofmass_1patt_vxy_angle_dif_fil(centerofmass_1patt_vxy_angle_dif_fil>angle_threshold_high) = nan;
        centerofmass_1patt_vxy_angle_dif_fil(centerofmass_1patt_vxy_angle_dif_fil<angle_threshold_low) = nan;
        
        % due to the ambiguity of phase values (and phase vector angles) 
        % close to the spiral centres, we assume 3x3 voxels centred at 
        % the spiral centre have a angle difference of 90 degrees, which
        % means they would always be included in the full-sized spiral
        
        centerofmass_1patt_round = round(centerofmass_1patt);
        for icol = centerofmass_1patt_round(1)-1:centerofmass_1patt_round(1)+1;
            for irow = centerofmass_1patt_round(2)-1:centerofmass_1patt_round(2)+1;
                if irow > 2 && icol > 2 && irow < size(temp1_phase,1)-1 && icol < size(temp1_phase,2)-1 % make sure no edges to cause error
                centerofmass_1patt_vxy_angle_dif_fil(irow,icol) = 0.5*pi; 
                end
            end
        end

        % incrementally increase the radias from the spiral centre, until 
        % the total count of voxels excluded from the full-sized spiral 
        % exceeds a threshold 
        % threshold = expansion radius * expansion threshold
        
        % exapansion threshold: stop increment expansion of spiral radius if exceeded
        expansion_threshold = 1;   
        for d = 2:0.5:100    % gradually growing expansion radius, with increments of 0.5 
            
            % radius filter, remove voxels with distances larger than expansion radius    
            centerofmass_1patt_vxy_abs_fil = centerofmass_1patt_vxy_abs; 
            centerofmass_1patt_vxy_abs_fil(centerofmass_1patt_vxy_abs_fil>d) = nan; 
            temp1 = centerofmass_1patt_vxy_abs_fil./centerofmass_1patt_vxy_abs_fil;
            temp1_count_notnan = nansum(temp1(:)); % count voxels left 
            
            % combine radius filters and angle difference filter
            temp2 = centerofmass_1patt_vxy_abs_fil .* centerofmass_1patt_vxy_angle_dif_fil; 
            temp2_count_notnan = nansum(temp2(:)./temp2(:)); % count voxels left

            % exapansion threshold: if the total number of voxels excluded 
            % exceeds the threshold (expansion radius * expansion threshold),
            % stop the expansion
            
            if temp1_count_notnan - temp2_count_notnan > expansion_threshold.*d; 
                temp2_01 = temp2./temp2;
                temp2_01(isnan(temp2_01)) = 0;
                if centerofmass_1patt_round(1)>size(temp1_phase,2)  % remove spiral at the edge of the cortical map
                    break
                end
                if centerofmass_1patt_round(2)>size(temp1_phase,1)  % remove spiral at the edge of the cortical map
                    break
                end
                if curlz(centerofmass_1patt_round(2),centerofmass_1patt_round(1),time) > 0
                spiral_filt_pos{ipatt,time} = sparse(temp2_01); %  anticlock-wise spiral, mark it as +1  
                end
                break
            end

        end
      end
end
% Detect candidate brain spirals via curl value calcualted from phase vector field: clockwise spiral only  
% clockwise and anticlockwise spirals are detected treated seperately
disp(['detecting full-sized clockwise brain spirals via phase vector field...'])    
curlz = [];
for time = 1:size(temp1_phase,3)      
    temp1_vx = Vx_flowmap_norm_phase(:,:,time);
    temp1_vy = Vy_flowmap_norm_phase(:,:,time);
    [x,y] = meshgrid(1:size(temp1_vx,2),1:size(temp1_vx,1));
    [curlz(:,:,time) cav(:,:,time)]  = curl(x,y,temp1_vx,temp1_vy);  % curl value of phase vector field    
end      
cav_filt_pos = curlz;
cav_filt_pos(cav_filt_pos>-1) = 0;     % detect only clockwise vorties
cav_filt_pos(isnan(cav_filt_pos)) = 0; % filter out voxels points with curl value > -1

% spiral core detection based on filtered curl values (requires function
% pattDetection_V4.m)

params.minPattTime = min_duration; % minimum duration of spiral (10 time steps)
params.minPattSize = 3;  % minimum size of spiral (3x3 matrix)

[WCentroids,absoluteTime,instantTotalPower,pattSize,patternIdx] = pattDetection_v4(cav_filt_pos,cav_filt_pos,params,0,'CenterOfMass_Amp');
clearvars instantTotalPower patternIdx pattSize temp1_vx temp1_vy   
 
% find the angle difference between center-originating vector (center of mass point to actual vector location) and phase vector

for ipatt = 1:size(absoluteTime,2) % select 1 spiral pattern detected   
      for time_withinpatt = 1:size(absoluteTime{ipatt},2);  % select 1 time point within the duration of the detect spiral pattern
        clearvars i
        y_tempalte = [];
        x_tempalte = [];
        centerofmass_1patt_vxy_angle_dif = [];
        centerofmass_1patt_vxy_angle_dif_fil = [];

        temp1_time = absoluteTime{ipatt};
        time = temp1_time(time_withinpatt); 
        
        % find the angles of all the phase vectors within the defined 
        % vector at the defined time 
        temp1_vx = Vx_flowmap_norm_phase(:,:,time); 
        temp1_vy = Vy_flowmap_norm_phase(:,:,time);
        temp1_vxy_angle = angle(temp1_vx + 1i.*temp1_vy); 
        
        % centre of mass position of the defined vector at the definded
        % time
        temp1_centerofmass_1patt = WCentroids{ipatt};
        centerofmass_1patt = temp1_centerofmass_1patt(time_withinpatt,:);
        
        % calculate the angles of spiral-centre-originated vectors toward all 
        % voxel positions across the flattened cortex
        for irow = 1:size(temp1_phase,1) 
            for icol = 1:size(temp1_phase,2) 
                y_tempalte(irow,:) = irow;
                x_tempalte(:,icol) = icol;   
            end
        end

        centerofmass_1patt_vx = round(x_tempalte - centerofmass_1patt(1));
        centerofmass_1patt_vy = round(y_tempalte - centerofmass_1patt(2));
        centerofmass_1patt_vxy_angle = angle(centerofmass_1patt_vx + 1i.*centerofmass_1patt_vy);
        
        % lengths of spiral-centre-originated vectors/ distances
        centerofmass_1patt_vxy_abs = sqrt(centerofmass_1patt_vx.*centerofmass_1patt_vx+centerofmass_1patt_vy.*centerofmass_1patt_vy);
        centerofmass_1patt_vxy_abs(centerofmass_1patt_vxy_abs==0) = 1;
        
        % angle differences between ortex-centre-originated vectors and
        % phase vectors at the same voxel
        centerofmass_1patt_vxy_angle_dif = centerofmass_1patt_vxy_angle - temp1_vxy_angle;

        for irow = 1:size(temp1_phase,1) 
            for icol = 1:size(temp1_phase,2) 
                temp1 = centerofmass_1patt_vxy_angle_dif(irow,icol);
                if temp1> pi        % limit the angle differences between -pi and pi
                    centerofmass_1patt_vxy_angle_dif(irow,icol) = centerofmass_1patt_vxy_angle_dif(irow,icol) - 2*pi;
                elseif temp1< -pi   % limit the angle differences between -pi and pi
                    centerofmass_1patt_vxy_angle_dif(irow,icol) = centerofmass_1patt_vxy_angle_dif(irow,icol) + 2*pi;
                end
            end
        end
        centerofmass_1patt_vxy_angle_dif = abs(centerofmass_1patt_vxy_angle_dif); 

        % angle difference threshold: phase vectors with angle differences 
        % >135 or <45 degrees from the centre-orginated-vector are excluded 
        % from the spiral, ideally angle difference should be 90 degrees
        % for a perfect spiral
        
        angle_threshold_low = pi./2 - 0.7854; %***********%0.7854; %45degrees- %0.5236; % 30degrees below %0.2618; % 15 degrees below
        angle_threshold_high = pi./2 + 0.7854; %*********%0.7854; %45degrees+ %0.2536; % 30degrees above %0.2618; % 15 degrees above
        centerofmass_1patt_vxy_angle_dif_fil = centerofmass_1patt_vxy_angle_dif;
        centerofmass_1patt_vxy_angle_dif_fil(centerofmass_1patt_vxy_angle_dif_fil>angle_threshold_high) = nan;
        centerofmass_1patt_vxy_angle_dif_fil(centerofmass_1patt_vxy_angle_dif_fil<angle_threshold_low) = nan;
        
        % due to the ambiguity of phase values (and phase vector angles) 
        % close to the spiral centres, we assume 3x3 voxels centred at 
        % the spiral centre have a angle difference of 90 degrees, which
        % means they would always be included in the full-sized spiral
        
        centerofmass_1patt_round = round(centerofmass_1patt);
        for icol = centerofmass_1patt_round(1)-1:centerofmass_1patt_round(1)+1;
            for irow = centerofmass_1patt_round(2)-1:centerofmass_1patt_round(2)+1;
                if irow > 2 && icol > 2 && irow < size(temp1_phase,1)-1 && icol < size(temp1_phase,2)-1 % make sure no edges to cause error
                centerofmass_1patt_vxy_angle_dif_fil(irow,icol) = 0.5*pi; 
                end
            end
        end

        % incrementally increase the radias from the spiral centre, until 
        % the total count of voxels excluded from the full-sized spiral 
        % exceeds a threshold 
        % threshold = expansion radius * expansion threshold
        
        % exapansion threshold: stop increment expansion of spiral radius if exceeded
        expansion_threshold = 1;   
        for d = 2:0.5:100    % gradually growing expansion radius, with increments of 0.5 
            
            % radius filter, remove voxels with distances larger than expansion radius    
            centerofmass_1patt_vxy_abs_fil = centerofmass_1patt_vxy_abs; 
            centerofmass_1patt_vxy_abs_fil(centerofmass_1patt_vxy_abs_fil>d) = nan; 
            temp1 = centerofmass_1patt_vxy_abs_fil./centerofmass_1patt_vxy_abs_fil;
            temp1_count_notnan = nansum(temp1(:)); % count voxels left 
            
            % combine radius filters and angle difference filter
            temp2 = centerofmass_1patt_vxy_abs_fil .* centerofmass_1patt_vxy_angle_dif_fil; 
            temp2_count_notnan = nansum(temp2(:)./temp2(:)); % count voxels left

            % exapansion threshold: if the total number of voxels excluded 
            % exceeds the threshold (expansion radius * expansion threshold),
            % stop the expansion
            
            if temp1_count_notnan - temp2_count_notnan > expansion_threshold.*d; 
                temp2_01 = temp2./temp2;
                temp2_01(isnan(temp2_01)) = 0;
                if centerofmass_1patt_round(1)>size(temp1_phase,2)  % remove spiral at the edge of the cortical map
                    break
                end
                if centerofmass_1patt_round(2)>size(temp1_phase,1) % remove spiral at the edge of the cortical map
                    break
                end
                if curlz(centerofmass_1patt_round(2),centerofmass_1patt_round(1),time) < 0
                spiral_filt_nega{ipatt,time} = sparse(-1.*temp2_01);  % if clock-wise spiral, mark it as -1 
                end
                break
            end

        end

      end
end 

%% filter spiral by size
size_threshold  = 25; %25; % 5x5 to be minimum size
spiral_filt_pos_sizefilt = spiral_filt_pos;
spiral_filt_nega_sizefilt = spiral_filt_nega;
for ipatt = 1:size(spiral_filt_pos,1)
    temp1_size_1patt = [];
    for t = 1:size(spiral_filt_pos,2)
        temp1 = full(spiral_filt_pos{ipatt,t});
        temp1_size = nansum(temp1(:));
        if  temp1_size ~= 0
            temp1_size_1patt = [temp1_size_1patt; temp1_size];
        end
    end
    temp1_size_1patt_maxsize = nanmax(abs(temp1_size_1patt(:)));
    if temp1_size_1patt_maxsize < size_threshold
        for t = 1:size(spiral_filt_pos,2)
            spiral_filt_pos_sizefilt{ipatt,t} = [];
        end
    end
end

for ipatt = 1:size(spiral_filt_nega,1)
    temp1_size_1patt = [];
    for t = 1:size(spiral_filt_nega,2)
        temp1 = full(spiral_filt_nega{ipatt,t});
        temp1_size = nansum(temp1(:));
        if  temp1_size ~= 0
            temp1_size_1patt = [temp1_size_1patt; temp1_size];
        end
    end
    temp1_size_1patt_maxsize = nanmax(abs(temp1_size_1patt(:)));
    if temp1_size_1patt_maxsize < size_threshold
        for t = 1:size(spiral_filt_nega,2)
            spiral_filt_nega_sizefilt{ipatt,t} = [];
        end
    end
end

% combine spiral template 
 for time = 1:window_size
spiral_template_math_listen = nan(175,251,2);

% positve only spiral
    for ipatt = 1:size(spiral_filt_pos,1);
        temp1 = spiral_filt_pos_sizefilt{ipatt,time};
        if nansum(temp1(:)) == 0;
            continue
        else
            temp_full = full(temp1);
            [row col] = find(temp_full);

            for i = 1:size(row,1)
                spiral_template_math_listen(row(i),col(i),1) = temp_full(row(i),col(i)) ;
            end
        end
    end
% negative only spiral
    for ipatt = 1:size(spiral_filt_nega,1);
        temp1 = spiral_filt_nega_sizefilt{ipatt,time};
        if nansum(temp1(:)) == 0;
            continue
        else
            temp_full = full(temp1);
            [row col] = find(temp_full);

            for i = 1:size(row,1)
                spiral_template_math_listen(row(i),col(i),2) = temp_full(row(i),col(i)) ;
            end
        end
    end
    spiral_template_math_listen(spiral_template_math_listen==2) = 1;    
    spiral_template_posi_nega = reshape(nansum(spiral_template_math_listen,3),[175,251]);
    spiral_template_posi_nega(isnan(spiral_template_posi_nega)) = 0; 
    spiral_template_posi_nega_accu(:,:,time) =  spiral_template_posi_nega;
    
 end



  %%  visualization
cd(main_folder)
if hemisphere == 1  
load('parcellation_template.mat')
elseif hemisphere == 2
load('parcellation_template22_RightBrain_subject1-100.mat')
parcellation_template_right = parcellation_template22_RightBrain_100sub(1:175,:,50);
end

downsample_rate = 2;

spiral_template_posi_nega_accu_filt = spiral_template_posi_nega_accu;
spiral_template_posi_nega_accu_filt(isnan(spiral_template_posi_nega_accu_filt)) = 0;

foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Analysis'];
cd(foldername)
filename = ['task_specific_spiral_distribution.mat'];
load(filename) 

% trial-averaged spiral distribution from spatiotemporally bandpass
% filtered fMRI signal (to be compared with spirals from unfiltered fMRI signal)
spiral_template_math_listen = zeros(175,251,window_size,size(spiral_distribution_math_listen_1to20,2));
for subject = 1:size(spiral_distribution_math_listen_1to20,2);
    for ipatt = 1:size(spiral_distribution_math_listen_1to20,1)
        temp1 = full(spiral_distribution_math_listen_1to20{ipatt,subject});
        if nansum(temp1(:))~=0
           spiral_template_math_listen(:,:,:,subject) = spiral_template_math_listen(:,:,:,subject) + temp1;
        end
    end
end
spiral_template_math_listen_avg = nanmean(spiral_template_math_listen,4);


spiral_template_math_answer = zeros(175,251,window_size,size(spiral_distribution_math_answer_1to20,2));
for subject = 1:size(spiral_distribution_math_answer_1to20,2);
    for ipatt = 1:size(spiral_distribution_math_answer_1to20,1)
        temp1 = full(spiral_distribution_math_answer_1to20{ipatt,subject});
        if nansum(temp1(:))~=0
           spiral_template_math_answer(:,:,:,subject) = spiral_template_math_answer(:,:,:,subject) + temp1;
        end
    end
end
spiral_template_math_answer_avg = nanmean(spiral_template_math_answer,4);

spiral_template_story_listen = zeros(175,251,window_size,size(spiral_distribution_story_listen_1to20,2));
for subject = 1:size(spiral_distribution_story_listen_1to20,2);
    for ipatt = 1:size(spiral_distribution_story_listen_1to20,1)
        temp1 = full(spiral_distribution_story_listen_1to20{ipatt,subject});
        if nansum(temp1(:))~=0
           spiral_template_story_listen(:,:,:,subject) = spiral_template_story_listen(:,:,:,subject) + temp1;
        end
    end
end
spiral_template_story_listen_avg = nanmean(spiral_template_story_listen,4);

spiral_template_story_answer = zeros(175,251,window_size,size(spiral_distribution_story_answer_1to20,2));

for subject = 1:size(spiral_distribution_story_answer_1to20,2);
    for ipatt = 1:size(spiral_distribution_story_answer_1to20,1)
        temp1 = full(spiral_distribution_story_answer_1to20{ipatt,subject});
        if nansum(temp1(:))~=0
           spiral_template_story_answer(:,:,:,subject) = spiral_template_story_answer(:,:,:,subject) + temp1;
        end
    end
end
spiral_template_story_answer_avg = nanmean(spiral_template_story_answer,4);


figure()
for t = 8:8 % t = 8  for listening tasks; t = 9 for answering tasks 
     temp1 = spiral_template_posi_nega_accu(:,:,t);
     temp1_pos = temp1;
     temp1_pos(temp1_pos<0) = 0;
     temp1_nega = temp1;
     temp1_nega(temp1_nega>0) = 0;
     
     [row_pos, col_pos] = find(temp1_pos);
     [row_nega, col_nega] = find(temp1_nega);

     % math listen and answer
     if listen_or_answer == 1
        temp2 = spiral_template_math_listen_avg(:,:,t); 
     elseif listen_or_answer == 2
        temp2 = spiral_template_math_answer_avg(:,:,t);  
     end
     
     % story listen and answer
%      if listen_or_answer == 1
%      temp2 = spiral_template_story_listen_avg(:,:,t);
%      elseif listen_or_answer == 2
%      temp2 = spiral_template_story_answer_avg(:,:,t);
%      end     

if hemisphere == 1  
 temp2 = temp2.*parcellation_template(1:175,1:251)./parcellation_template(1:175,1:251);
elseif hemisphere == 2
 temp2 = temp2.*parcellation_template_right(1:175,1:251)./parcellation_template_right(1:175,1:251);
end
min = nanmean(temp2(:));
max = nanmax(temp2(:));
if min > max 
    min_adj = min;
    max_adj = -1.*min;
else
    min_adj = -1.*max;
    max_adj = max;
end     
    pcolor(temp2)
    shading interp
    caxis([min_adj,max_adj])
%         caxis([-0.15,0.15])
    xlim([1,251])
    ylim([1,175])
    colorbar 
    colormap jet
    hold on
     for parcellation_ID = 1:22
         if hemisphere == 1
            parcellation_template_1par = parcellation_template(1:175,1:251);
         elseif hemisphere == 2
            parcellation_template_1par = parcellation_template22_RightBrain_100sub(1:175,1:251,1);
         end
        parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
        parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;
        B = bwboundaries(parcellation_template_1par,'noholes');
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2),boundary(:,1),'-','linewidth',1,'color',[0.5,0.5,0.5])
        end
     end 
     hold off   
     hold on
     scatter(col_pos,row_pos,2,[0 0 0],'filled')
     scatter(col_nega,row_nega,2,[1 1 1],'filled')          
     set(gca,'ydir','normal')
     hold off
    title(['task-evoked unfiltered fMRI signal, math listen,t=',num2str(t)])

end

      sgtitle(['max size>',num2str(size_threshold),', 1s movmean of trial averaged amplitude, bandpass 83.62-7.45']) 



foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Analysis/'];
filename = ['task_evoked_unfiltered_fMRI_signal.mat'];
save([foldername,filename],'spiral_template_posi_nega_accu','row_pos','col_pos'...
    ,'col_nega','row_nega','spiral_template_math_listen_avg','spiral_template_math_answer_avg'...
    ,'spiral_template_story_listen_avg','spiral_template_story_answer_avg'...
    ,'data_math_listen_avg','data_math_answer_avg','temp1_phase'...
    ,'Vx_flowmap_norm_phase','Vy_flowmap_norm_phase')


%% visualziation: zoomed in view of the task-evoked unfiltered fMRI signal 


if listen_or_answer == 1
    data_math_listen_avg_min = nanmin(data_math_listen_avg,[],3);
    data_math_listen_avg_max = nanmax(data_math_listen_avg,[],3);
    data_math_listen_avg_minmaxnorm = (data_math_listen_avg - data_math_listen_avg_min)./(data_math_listen_avg_max - data_math_listen_avg_min);
elseif listen_or_answer == 2
    data_math_answer_avg_min = nanmin(data_math_answer_avg,[],3);
    data_math_answer_avg_max = nanmax(data_math_answer_avg,[],3);
    data_math_answer_avg_minmaxnorm = (data_math_answer_avg - data_math_answer_avg_min)./(data_math_answer_avg_max - data_math_answer_avg_min);
end
figure()
for t = 1:20
    subplot(4,5,t)
    if listen_or_answer == 1
    pcolor(data_math_listen_avg_minmaxnorm(:,:,t))
    elseif listen_or_answer == 2
    pcolor(data_math_answer_avg_minmaxnorm(:,:,t))
    end
    shading interp
    ylim([128,158]) % PCC-Left listen/answer
    xlim([175,205])
%     ylim([130,160]) % PCC-R answer/listen
%     xlim([65,95])  
% if hemisphere == 1
%     if listen_or_answer == 1
%     ylim([90,120]) % Motor-L-listen
%     xlim([110,140])   
%     elseif listen_or_answer == 2
%     ylim([95,125]) % Motor-L-listen
%     xlim([110,140])   
%     end
% elseif hemisphere == 2
% %     ylim([105,135]) % Motor-R-listen
%     ylim([95,125]) % Motor-R-listen and answer
%     xlim([120,150])   
% end
hold on
for parcellation_ID = 1:22
    if hemisphere == 1
    parcellation_template_1par = parcellation_template(1:175,1:251);
elseif hemisphere == 2
    parcellation_template_1par = parcellation_template22_RightBrain_100sub(1:175,1:251,100);
end
parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;
B = bwboundaries(parcellation_template_1par,'noholes');
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2),boundary(:,1),'-','linewidth',1,'color',[1,1,1])
end
end 
hold off   
% pause(0.1)
end

% Motor_L_listen_xylim = ['y 90 120 x 110 140 '];
% Motor_R_listen_xylim = ['y 95 125 x 120 150'];
% Motor_L_answer_xylim = ['y 90 125 x 110 140'];
% Motor_R_answer_xylim = ['y 90 125 x 110 150'];
% PCC_L_answer_xylim = ['y 128 158 x 175 205'];
% PCC_R_answer_xylim = ['y 130 160 x 65 95'];
% PCC_R_listen_xylim = ['y 130 160 x 65 95'];

%% phase vecotr field as streamlines (activity flows)

Vx_flowmap_norm_phase_movmean = movmean(Vx_flowmap_norm_phase,10,3);
Vy_flowmap_norm_phase_movmean = movmean(Vy_flowmap_norm_phase,10,3);

figure()
for t = 8:8
    temp1_x = Vx_flowmap_norm_phase_movmean(:,:,t);
    temp1_y = Vy_flowmap_norm_phase_movmean(:,:,t);
    temp1_x_norm = temp1_x./sqrt(temp1_x.^2 + temp1_y.^2);
    temp1_y_norm = temp1_y./sqrt(temp1_x.^2 + temp1_y.^2); 
    
    temp1 = spiral_template_posi_nega_accu_filt(:,:,t);

    imagesc(temp1)
    set(gca,'ydir','normal')
    colorbar
    
%     title(['Amp.Grad.,MathAns.,t=',num2str(t)])
%     title(['Amp.,MathAns.,t=',num2str(t)])
%     title(['Amp.,MathLis.,t=',num2str(t)])
%     title(['Amp.Grad.,MathLis.,t=',num2str(t)])
%     title(['Amp.,MathLis.,t=',num2str(t)])
%     title(['Amp.,StoryLis.,t=',num2str(t)])
%     title(['Trial averaged Amplitude, Math Listen, ',num2str(t.*0.72),'s after onset'])
    hold on
    for parcellation_ID = 1:22
        if hemisphere == 1
            parcellation_template_1par = parcellation_template;
        elseif hemisphere == 2
            parcellation_template_1par = parcellation_template22_RightBrain_100sub(:,:,100);
        end
    parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
    parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;
    B = bwboundaries(parcellation_template_1par,'noholes');
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2),boundary(:,1),'-','linewidth',2,'color',[0.5,0.5,0.5])
    end
    end 
    hold off   
    hold on
%     Q1 = quiver(temp1_x_norm(1:175,1:251),temp1_y_norm(1:175,1:251),1,'color',[0,0,0]);
%     Q1.LineWidth =2; 
    h = streamslice([1:251],[1:175],temp1_x,temp1_y,110,'k');
%     h = streamslice([1:251],[1:175],temp1_x,temp1_y,40,'k');
    % left brain ONLY
%     ylim([138,148]) % PCC-Left answer
%     xlim([185,195])
    ylim([137,147]) % PCC-Left listen
    xlim([183,193])
%     ylim([90,120])
%     xlim([110,140])    
%     ylim([90,120])
%     xlim([170,200])    
    % right brain ONLY
%     ylim([139,149]) %PCC-R_listen
%     xlim([75,85])    
%     ylim([139,149]) %PCC-R_answer
%     xlim([72,82])      
%     ylim([115,145])
%     xlim([60,90])   
%     ylim([98,128])
%     xlim([125,155])          
%     ylim([90,110])
%     xlim([75,95]) 
    hold off
end
% temp1_x_window = temp1_x(139:149,72:82); %PCC-R_answer
% temp1_y_window = temp1_y(139:149,72:82);
% temp1_x_window = temp1_x(138:148,185:195);  % PCC-Left answer
% temp1_y_window = temp1_y(138:148,185:195);
% temp1_x_window = temp1_x(137:147,183:193);  % PCC-Left listen
% temp1_y_window = temp1_y(137:147,183:193);
% temp1_x_window = temp1_x(139:149,75:85);  % PCC-R listen
% temp1_y_window = temp1_y(139:149,75:85);


end
    