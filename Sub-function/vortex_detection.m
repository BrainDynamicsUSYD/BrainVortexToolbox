function [vortex_filt_nega,vortex_filt_pos,subject] = vortex_detection(subject)
%% detect full-sized brain vortices from pre-processed fMRI data
disp(['initiating process...'])    
disp(['starting subject ',num2str(subject)])    

vortex_filt_nega = [];
vortex_filt_pos = [];

% load preprocessed fMRI data files
disp(['loading preprocessed fMRI data files...'])    

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
disp(['producing phase gradient velocity field baesd on phase map...'])    

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

%% Detect candidate brain vortices via curl value calcualted from phase vector field: anticlockwise vortex only  
% clockwise and anticlockwise vortices are detected treated seperately
disp(['detecting full-sized anticlockwise brain vortices via phase vector field...'])    
  
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
 
% find the angle difference between center-originating vector (center of mass point to actual vector location) and phase vector

load parcellation_template % load parcellation template to confirm the edge of flattened cortical map 
parcellation_template_01 = parcellation_template./parcellation_template;

for ipatt = 1:size(absoluteTime,2) % select 1 vortex pattern detected   
      for time_withinpatt = 1:size(absoluteTime{ipatt},2);  % select 1 time point within the duration of the detect vortex pattern
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
        temp1_vxy_angle = angle(temp1_vx + i.*temp1_vy); 
        
        % centre of mass position of the defined vector at the definded
        % time
        temp1_centerofmass_1patt = WCentroids{ipatt};
        centerofmass_1patt = temp1_centerofmass_1patt(time_withinpatt,:);

        
        % calculate the angles of vortex-centre-originated vectors toward all 
        % voxel positions across the flattened cortex
        for irow = 1:176
            for icol = 1:251
                y_tempalte(irow,:) = irow;
                x_tempalte(:,icol) = icol;   
            end
        end
        x_tempalte = x_tempalte.*parcellation_template_01;
        y_tempalte = y_tempalte.*parcellation_template_01;

        centerofmass_1patt_vx = round(x_tempalte - centerofmass_1patt(1));
        centerofmass_1patt_vy = round(y_tempalte - centerofmass_1patt(2));
        centerofmass_1patt_vxy_angle = angle(centerofmass_1patt_vx + i.*centerofmass_1patt_vy);
        
        % lengths of vortex-centre-originated vectors/ distances
        centerofmass_1patt_vxy_abs = sqrt(centerofmass_1patt_vx.*centerofmass_1patt_vx+centerofmass_1patt_vy.*centerofmass_1patt_vy);
        centerofmass_1patt_vxy_abs(centerofmass_1patt_vxy_abs==0) = 1;
        
        % angle differences between ortex-centre-originated vectors and
        % phase vectors at the same voxel
        centerofmass_1patt_vxy_angle_dif = centerofmass_1patt_vxy_angle - temp1_vxy_angle;
        
        
        for irow = 1:176
            for icol = 1:251
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
        % from the vortex, ideally angle difference should be 90 degrees
        % for a perfect vortex
        
        angle_threshold_low = pi./2 - 0.7854; %***********%0.7854; %45degrees- %0.5236; % 30degrees below %0.2618; % 15 degrees below
        angle_threshold_high = pi./2 + 0.7854; %*********%0.7854; %45degrees+ %0.2536; % 30degrees above %0.2618; % 15 degrees above
        centerofmass_1patt_vxy_angle_dif_fil = centerofmass_1patt_vxy_angle_dif;
        centerofmass_1patt_vxy_angle_dif_fil(centerofmass_1patt_vxy_angle_dif_fil>angle_threshold_high) = nan;
        centerofmass_1patt_vxy_angle_dif_fil(centerofmass_1patt_vxy_angle_dif_fil<angle_threshold_low) = nan;
        
        % due to the ambiguity of phase values (and phase vector angles) 
        % close to the vortex centres, we assume 3x3 voxels centred at 
        % the vortex centre have a angle difference of 90 degrees, which
        % means they would always be included in the full-sized vortex
        
        centerofmass_1patt_round = round(centerofmass_1patt);
        for icol = centerofmass_1patt_round(1)-1:centerofmass_1patt_round(1)+1;
            for irow = centerofmass_1patt_round(2)-1:centerofmass_1patt_round(2)+1;
                centerofmass_1patt_vxy_angle_dif_fil(irow,icol) = 0.5*pi; 
            end
        end

        % incrementally increase the radias from the vortex centre, until 
        % the total count of voxels excluded from the full-sized vortex 
        % exceeds a threshold 
        % threshold = expansion radius * expansion threshold
        
        % exapansion threshold: stop increment expansion of vortex radius if exceeded
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
                if centerofmass_1patt_round(1)>251 % remove vortex at the edge of the cortical map
                    break
                end
                if centerofmass_1patt_round(2)>176 % remove vortex at the edge of the cortical map
                    break
                end
                if curlz(centerofmass_1patt_round(2),centerofmass_1patt_round(1),time) > 0
                vortex_filt_pos{ipatt,time} = sparse(temp2_01); %  anticlock-wise vortex, mark it as +1  
                end
                break
            end

        end
      end
end
%% Detect candidate brain vortices via curl value calcualted from phase vector field: clockwise vortex only  
% clockwise and anticlockwise vortices are detected treated seperately
disp(['detecting full-sized clockwise brain vortices via phase vector field...'])    
  
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
 
% find the angle difference between center-originating vector (center of mass point to actual vector location) and phase vector

load parcellation_template % load parcellation template to confirm the edge of flattened cortical map 
parcellation_template_01 = parcellation_template./parcellation_template;

for ipatt = 1:size(absoluteTime,2) % select 1 vortex pattern detected   
      for time_withinpatt = 1:size(absoluteTime{ipatt},2);  % select 1 time point within the duration of the detect vortex pattern
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
        temp1_vxy_angle = angle(temp1_vx + i.*temp1_vy); 
        
        % centre of mass position of the defined vector at the definded
        % time
        temp1_centerofmass_1patt = WCentroids{ipatt};
        centerofmass_1patt = temp1_centerofmass_1patt(time_withinpatt,:);

        
        % calculate the angles of vortex-centre-originated vectors toward all 
        % voxel positions across the flattened cortex
        for irow = 1:176
            for icol = 1:251
                y_tempalte(irow,:) = irow;
                x_tempalte(:,icol) = icol;   
            end
        end
        x_tempalte = x_tempalte.*parcellation_template_01;
        y_tempalte = y_tempalte.*parcellation_template_01;

        centerofmass_1patt_vx = round(x_tempalte - centerofmass_1patt(1));
        centerofmass_1patt_vy = round(y_tempalte - centerofmass_1patt(2));
        centerofmass_1patt_vxy_angle = angle(centerofmass_1patt_vx + i.*centerofmass_1patt_vy);
        
        % lengths of vortex-centre-originated vectors/ distances
        centerofmass_1patt_vxy_abs = sqrt(centerofmass_1patt_vx.*centerofmass_1patt_vx+centerofmass_1patt_vy.*centerofmass_1patt_vy);
        centerofmass_1patt_vxy_abs(centerofmass_1patt_vxy_abs==0) = 1;
        
        % angle differences between ortex-centre-originated vectors and
        % phase vectors at the same voxel
        centerofmass_1patt_vxy_angle_dif = centerofmass_1patt_vxy_angle - temp1_vxy_angle;
        
        
        for irow = 1:176
            for icol = 1:251
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
        % from the vortex, ideally angle difference should be 90 degrees
        % for a perfect vortex
        
        angle_threshold_low = pi./2 - 0.7854; %***********%0.7854; %45degrees- %0.5236; % 30degrees below %0.2618; % 15 degrees below
        angle_threshold_high = pi./2 + 0.7854; %*********%0.7854; %45degrees+ %0.2536; % 30degrees above %0.2618; % 15 degrees above
        centerofmass_1patt_vxy_angle_dif_fil = centerofmass_1patt_vxy_angle_dif;
        centerofmass_1patt_vxy_angle_dif_fil(centerofmass_1patt_vxy_angle_dif_fil>angle_threshold_high) = nan;
        centerofmass_1patt_vxy_angle_dif_fil(centerofmass_1patt_vxy_angle_dif_fil<angle_threshold_low) = nan;
        
        % due to the ambiguity of phase values (and phase vector angles) 
        % close to the vortex centres, we assume 3x3 voxels centred at 
        % the vortex centre have a angle difference of 90 degrees, which
        % means they would always be included in the full-sized vortex
        
        centerofmass_1patt_round = round(centerofmass_1patt);
        for icol = centerofmass_1patt_round(1)-1:centerofmass_1patt_round(1)+1;
            for irow = centerofmass_1patt_round(2)-1:centerofmass_1patt_round(2)+1;
                centerofmass_1patt_vxy_angle_dif_fil(irow,icol) = 0.5*pi; 
            end
        end

        % incrementally increase the radias from the vortex centre, until 
        % the total count of voxels excluded from the full-sized vortex 
        % exceeds a threshold 
        % threshold = expansion radius * expansion threshold
        
        % exapansion threshold: stop increment expansion of vortex radius if exceeded
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
                if centerofmass_1patt_round(1)>251 % remove vortex at the edge of the cortical map
                    break
                end
                if centerofmass_1patt_round(2)>176 % remove vortex at the edge of the cortical map
                    break
                end
                if curlz(centerofmass_1patt_round(2),centerofmass_1patt_round(1),time) < 0
                vortex_filt_nega{ipatt,time} = sparse(-1.*temp2_01);  % if clock-wise vortex, mark it as -1 
                end
                break
            end

        end
      end
end
        save_folder = [pwd,'/Sample Data/'];
        save([save_folder,'vortex_bandpass_29.35_14.93_subject_',num2str(subject),'.mat'],'vortex_filt_nega','vortex_filt_pos','subject')        
        disp(['finishing subject ',num2str(subject)])    
    
    
end
