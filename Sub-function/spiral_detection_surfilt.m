function     [temp1_compatibility_ratio_pos_nega_real_avg,temp1_compatibility_ratio_pos_nega_sur_avg] = spiral_detection_surfilt(subject,DataIn_smooth,DataIn_unsmooth, DataIn_sur_smooth, DataIn_sur_unsmooth,main_folder,flagRest,flagSmooth,flagTask,hemisphere)

cd(main_folder)
expansion_threshold = 1;
phase_df_threshold = pi/6;
%% surrogate data
% ensure surrogate data (raw & badpass) is compatible with the flattened
% cortical map
if hemisphere == 1
    load('parcellation_template.mat')
elseif hemisphere == 2
    load('parcellation_template22_RightBrain_subject1-100.mat')
    parcellation_template = parcellation_template22_RightBrain_100sub(:,:,1);
end

DataIn_sur_smooth = DataIn_sur_smooth.*parcellation_template(1:175,1:251)./parcellation_template(1:175,1:251);
DataIn_sur_unsmooth = DataIn_sur_unsmooth.*parcellation_template(1:175,1:251)./parcellation_template(1:175,1:251);

temp1_smooth_phase_sur = nan(size(DataIn_sur_smooth));
temp1_raw_phase_sur = nan(size(DataIn_sur_unsmooth));
for irow = 1:size(temp1_smooth_phase_sur,1)
    for icol = 1:size(temp1_smooth_phase_sur,2)
        temp1 = DataIn_sur_smooth(irow,icol,:);
        if nansum(temp1(:))~=0
        temp1_smooth_phase_sur(irow,icol,:) = angle(hilbert(temp1(:)));
        end
        temp1 = DataIn_sur_unsmooth(irow,icol,:);
        if nansum(temp1(:))~=0
        temp1_raw_phase_sur(irow,icol,:) = angle(hilbert(temp1(:)));
        end        
        
    end
end

% phase gradient (vector) field

phaseSig = temp1_smooth_phase_sur;
spiral_filt_nega_sur = []; 
spiral_filt_pos_sur = [];
Vx_flowmap_norm_phase = [];
Vy_flowmap_norm_phase = [];
spiral_size_nega_sur= [];
spiral_size_pos_sur = [];
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

% anticlockwise spirals
% Detect candidate brain spirals via curl value calcualted from phase vector field: anticlockwise spiral only  
% clockwise and anticlockwise spirals are detected treated seperately
curlz = [];
cav = [];  
for time = 1:size(Vx_flowmap_norm_phase,3)     
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

params.minPattTime = 1; % minimum duration of spiral (10 time steps)
params.minPattSize = 3;  % minimum size of spiral (3x3 matrix)

[WCentroids,absoluteTime,instantTotalPower,pattSize,patternIdx] = pattDetection_v4(cav_filt_pos,cav_filt_pos,params,0,'CenterOfMass_Amp');
clearvars instantTotalPower patternIdx pattSize temp1_vx temp1_vy   
 
% find the centre ONLY positions in x/y coordinates
for ipatt = 1:size(WCentroids,2)
    temp1_absoluteTime = absoluteTime{ipatt};
    temp1_xy = WCentroids{ipatt};
    temp1_x = round(temp1_xy(:,1),2);
    temp1_y = round(temp1_xy(:,2),2);
    for time = 1:size(temp1_absoluteTime,2)
        spiral_filt_pos_centreONLY_sur{ipatt,temp1_absoluteTime(time)} = [temp1_x(time) temp1_y(time)];
    end
end
        
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
        temp1_vxy_angle = angle(temp1_vx + i.*temp1_vy); 
        
        % centre of mass position of the defined vector at the definded
        % time
        temp1_centerofmass_1patt = WCentroids{ipatt};
        centerofmass_1patt = temp1_centerofmass_1patt(time_withinpatt,:);
        
        % calculate the angles of spiral-centre-originated vectors toward all 
        % voxel positions across the flattened cortex
        for irow = 1:size(Vx_flowmap_norm_phase,1) 
            for icol = 1:size(Vx_flowmap_norm_phase,2) 
                y_tempalte(irow,:) = irow;
                x_tempalte(:,icol) = icol;   
            end
        end
        centerofmass_1patt_vx = round(x_tempalte - centerofmass_1patt(1));
        centerofmass_1patt_vy = round(y_tempalte - centerofmass_1patt(2));
        centerofmass_1patt_vxy_angle = angle(centerofmass_1patt_vx + i.*centerofmass_1patt_vy);
        
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
        angle_threshold_low = pi./2 - 0.7854; %***********%45degrees
        angle_threshold_high = pi./2 + 0.7854; %*********%45degrees
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
                if irow > 2 && icol > 2 && irow < size(Vx_flowmap_norm_phase,1)-1 && icol < size(Vx_flowmap_norm_phase,2)-1 % make sure no edges to cause error
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
                if centerofmass_1patt_round(1)>size(Vx_flowmap_norm_phase,2)  % remove spiral at the edge of the cortical map
                    break
                end
                if centerofmass_1patt_round(2)>size(Vx_flowmap_norm_phase,1)  % remove spiral at the edge of the cortical map
                    break
                end
                if curlz(centerofmass_1patt_round(2),centerofmass_1patt_round(1),time) > 0
                spiral_filt_pos_sur{ipatt,time} = sparse(temp2_01); %  anticlock-wise spiral, mark it as +1  
                end
                break
            end

        end
        spiral_size_pos_sur{ipatt,time} = d; % final size of each spiral at each time
      end
end

% clockwise spirals (repeat the same process as with anticlockwise spirals)
curlz = [];
cav = [];
for time = 1:size(Vx_flowmap_norm_phase,3)      
    temp1_vx = Vx_flowmap_norm_phase(:,:,time);
    temp1_vy = Vy_flowmap_norm_phase(:,:,time);
    [x,y] = meshgrid(1:size(temp1_vx,2),1:size(temp1_vx,1));
    [curlz(:,:,time) cav(:,:,time)]  = curl(x,y,temp1_vx,temp1_vy);  % curl value of phase vector field    
end      
cav_filt_nega = curlz; % detect only clockwise vorties
cav_filt_nega(cav_filt_nega>-1) = 0;     
cav_filt_nega(isnan(cav_filt_nega)) = 0; % filter out voxels points with curl value > -1

% spiral core detection based on filtered curl values (requires function
% pattDetection_V4.m)

params.minPattTime = 1; % minimum duration of spiral (10 time steps)
params.minPattSize = 3;  % minimum size of spiral (3x3 matrix)

[WCentroids,absoluteTime,instantTotalPower,pattSize,patternIdx] = pattDetection_v4(cav_filt_nega,cav_filt_nega,params,0,'CenterOfMass_Amp');
clearvars instantTotalPower patternIdx pattSize temp1_vx temp1_vy   
 
% find the centre ONLY positions in x/y coordinates
for ipatt = 1:size(WCentroids,2)
    temp1_absoluteTime = absoluteTime{ipatt};
    temp1_xy = WCentroids{ipatt};
    temp1_x = round(temp1_xy(:,1),2);
    temp1_y = round(temp1_xy(:,2),2);
    for time = 1:size(temp1_absoluteTime,2)
        spiral_filt_nega_centreONLY_sur{ipatt,temp1_absoluteTime(time)} = [temp1_x(time) temp1_y(time)];
    end
end

% find the angle difference between center-originating vector (center of mass point to actual vector location) and phase vector

% load parcellation_template % load parcellation template to confirm the edge of flattened cortical map 
% parcellation_template_01 = parcellation_template./parcellation_template;


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
        temp1_vxy_angle = angle(temp1_vx + i.*temp1_vy); 
        
        % centre of mass position of the defined vector at the definded
        % time
        temp1_centerofmass_1patt = WCentroids{ipatt};
        centerofmass_1patt = temp1_centerofmass_1patt(time_withinpatt,:);
        
        
        % calculate the angles of spiral-centre-originated vectors toward all 
        % voxel positions across the flattened cortex
        for irow = 1:size(Vx_flowmap_norm_phase,1) 
            for icol = 1:size(Vx_flowmap_norm_phase,2) 
                y_tempalte(irow,:) = irow;
                x_tempalte(:,icol) = icol;   
            end
        end
        centerofmass_1patt_vx = round(x_tempalte - centerofmass_1patt(1));
        centerofmass_1patt_vy = round(y_tempalte - centerofmass_1patt(2));
        centerofmass_1patt_vxy_angle = angle(centerofmass_1patt_vx + i.*centerofmass_1patt_vy);
        
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
        
        angle_threshold_low = pi./2 - 0.7854; %***********%45degrees 
        angle_threshold_high = pi./2 + 0.7854; %*********%45degrees
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
                if irow > 2 && icol > 2 && irow < size(Vx_flowmap_norm_phase,1)-1 && icol < size(Vx_flowmap_norm_phase,2)-1 % make sure no edges to cause error
                centerofmass_1patt_vxy_angle_dif_fil(irow,icol) = 0.5*pi; 
                end
            end
        end

        % incrementally increase the radias from the spiral centre, until 
        % the total count of voxels excluded from the full-sized spiral 
        % exceeds a threshold 
        % threshold = expansion radius * expansion threshold
        
        % exapansion threshold: stop increment expansion of spiral radius if exceeded
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
                if centerofmass_1patt_round(1)>size(Vx_flowmap_norm_phase,2)  % remove spiral at the edge of the cortical map
                    break
                end
                if centerofmass_1patt_round(2)>size(Vx_flowmap_norm_phase,1) % remove spiral at the edge of the cortical map
                    break
                end
                if curlz(centerofmass_1patt_round(2),centerofmass_1patt_round(1),time) < 0
                spiral_filt_nega_sur{ipatt,time} = sparse(-1.*temp2_01);  % if clock-wise spiral, mark it as -1 
                end
                break
            end
        end
                spiral_size_nega_sur{ipatt,time} = d; % final size of each spiral at each time
      end
end


%% real data



temp1_smooth_phase_real = nan(size(DataIn_smooth(1:175,:,:)));
temp1_raw_phase_real = nan(size(DataIn_unsmooth(1:175,:,:)));
for irow = 1:size(temp1_smooth_phase_real,1)
    for icol = 1:size(temp1_smooth_phase_real,2)
        temp1 = DataIn_smooth(irow,icol,:);
        if nansum(temp1(:))~=0
        temp1_smooth_phase_real(irow,icol,:) = angle(hilbert(temp1(:)));
        end
        temp1 = DataIn_unsmooth(irow,icol,:);
        if nansum(temp1(:))~=0
        temp1_raw_phase_real(irow,icol,:) = angle(hilbert(temp1(:)));
        end        
    end
end

% phase gradient (vector) field

phaseSig = temp1_smooth_phase_real;
spiral_filt_nega_real = [];
spiral_filt_pos_real = [];
Vx_flowmap_norm_phase = [];
Vy_flowmap_norm_phase = [];
spiral_size_nega_real= [];
spiral_size_pos_real = [];
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

% anticlockwise sprial
% Detect candidate brain spirals via curl value calcualted from phase vector field: anticlockwise spiral only  
% clockwise and anticlockwise spirals are detected treated seperately
curlz = [];
cav = [];
for time = 1:size(Vx_flowmap_norm_phase,3)     
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

params.minPattTime = 1; % minimum duration of spiral (10 time steps)
params.minPattSize = 3;  % minimum size of spiral (3x3 matrix)

[WCentroids,absoluteTime,instantTotalPower,pattSize,patternIdx] = pattDetection_v4(cav_filt_pos,cav_filt_pos,params,0,'CenterOfMass_Amp');
clearvars instantTotalPower patternIdx pattSize temp1_vx temp1_vy   
 
% find the centre ONLY positions in x/y coordinates
for ipatt = 1:size(WCentroids,2)
    temp1_absoluteTime = absoluteTime{ipatt};
    temp1_xy = WCentroids{ipatt};
    temp1_x = round(temp1_xy(:,1),2);
    temp1_y = round(temp1_xy(:,2),2);
    for time = 1:size(temp1_absoluteTime,2)
        spiral_filt_pos_centreONLY_real{ipatt,temp1_absoluteTime(time)} = [temp1_x(time) temp1_y(time)];
    end
end
        
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
        temp1_vxy_angle = angle(temp1_vx + i.*temp1_vy); 
        
        % centre of mass position of the defined vector at the definded
        % time
        temp1_centerofmass_1patt = WCentroids{ipatt};
        centerofmass_1patt = temp1_centerofmass_1patt(time_withinpatt,:);
        
        % calculate the angles of spiral-centre-originated vectors toward all 
        % voxel positions across the flattened cortex
        for irow = 1:size(Vx_flowmap_norm_phase,1) 
            for icol = 1:size(Vx_flowmap_norm_phase,2) 
                y_tempalte(irow,:) = irow;
                x_tempalte(:,icol) = icol;   
            end
        end
        centerofmass_1patt_vx = round(x_tempalte - centerofmass_1patt(1));
        centerofmass_1patt_vy = round(y_tempalte - centerofmass_1patt(2));
        centerofmass_1patt_vxy_angle = angle(centerofmass_1patt_vx + i.*centerofmass_1patt_vy);
        
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
                if irow > 2 && icol > 2 && irow < size(Vx_flowmap_norm_phase,1)-1 && icol < size(Vx_flowmap_norm_phase,2)-1 % make sure no edges to cause error
                centerofmass_1patt_vxy_angle_dif_fil(irow,icol) = 0.5*pi; 
                end
            end
        end

        % incrementally increase the radias from the spiral centre, until 
        % the total count of voxels excluded from the full-sized spiral 
        % exceeds a threshold 
        % threshold = expansion radius * expansion threshold
        
        % exapansion threshold: stop increment expansion of spiral radius if exceeded
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
                if centerofmass_1patt_round(1)>size(Vx_flowmap_norm_phase,2)  % remove spiral at the edge of the cortical map
                    break
                end
                if centerofmass_1patt_round(2)>size(Vx_flowmap_norm_phase,1)  % remove spiral at the edge of the cortical map
                    break
                end
                if curlz(centerofmass_1patt_round(2),centerofmass_1patt_round(1),time) > 0
                spiral_filt_pos_real{ipatt,time} = sparse(temp2_01); %  anticlock-wise spiral, mark it as +1  
                end
                break
            end

        end
        spiral_size_pos_real{ipatt,time} = d; % final size of each spiral at each time
      end
end

% clockwise spirals
curlz = [];
cav = [];  
for time = 1:size(Vx_flowmap_norm_phase,3)      
    temp1_vx = Vx_flowmap_norm_phase(:,:,time);
    temp1_vy = Vy_flowmap_norm_phase(:,:,time);
    [x,y] = meshgrid(1:size(temp1_vx,2),1:size(temp1_vx,1));
    [curlz(:,:,time) cav(:,:,time)]  = curl(x,y,temp1_vx,temp1_vy);  % curl value of phase vector field    
end      
cav_filt_nega = curlz;
cav_filt_nega(cav_filt_nega>-1) = 0;     % detect only clockwise vorties
cav_filt_nega(isnan(cav_filt_nega)) = 0; % filter out voxels points with curl value > -1

% spiral core detection based on filtered curl values (requires function
% pattDetection_V4.m)

params.minPattTime = 1; % minimum duration of spiral (10 time steps)
params.minPattSize = 3;  % minimum size of spiral (3x3 matrix)

[WCentroids,absoluteTime,instantTotalPower,pattSize,patternIdx] = pattDetection_v4(cav_filt_nega,cav_filt_nega,params,0,'CenterOfMass_Amp');
clearvars instantTotalPower patternIdx pattSize temp1_vx temp1_vy   
 
% find the centre ONLY positions in x/y coordinates
for ipatt = 1:size(WCentroids,2)
    temp1_absoluteTime = absoluteTime{ipatt};
    temp1_xy = WCentroids{ipatt};
    temp1_x = round(temp1_xy(:,1),2);
    temp1_y = round(temp1_xy(:,2),2);
    for time = 1:size(temp1_absoluteTime,2)
        spiral_filt_nega_centreONLY_real{ipatt,temp1_absoluteTime(time)} = [temp1_x(time) temp1_y(time)];
    end
end

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
        temp1_vxy_angle = angle(temp1_vx + i.*temp1_vy); 
        
        % centre of mass position of the defined vector at the definded
        % time
        temp1_centerofmass_1patt = WCentroids{ipatt};
        centerofmass_1patt = temp1_centerofmass_1patt(time_withinpatt,:);
        
        % calculate the angles of spiral-centre-originated vectors toward all 
        % voxel positions across the flattened cortex
        for irow = 1:size(Vx_flowmap_norm_phase,1) 
            for icol = 1:size(Vx_flowmap_norm_phase,2) 
                y_tempalte(irow,:) = irow;
                x_tempalte(:,icol) = icol;   
            end
        end
        centerofmass_1patt_vx = round(x_tempalte - centerofmass_1patt(1));
        centerofmass_1patt_vy = round(y_tempalte - centerofmass_1patt(2));
        centerofmass_1patt_vxy_angle = angle(centerofmass_1patt_vx + i.*centerofmass_1patt_vy);
        
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
        
        angle_threshold_low = pi./2 - 0.7854; %***********%45degrees
        angle_threshold_high = pi./2 + 0.7854; %*********%45degrees
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
                if irow > 2 && icol > 2 && irow < size(Vx_flowmap_norm_phase,1)-1 && icol < size(Vx_flowmap_norm_phase,2)-1 % make sure no edges to cause error
                centerofmass_1patt_vxy_angle_dif_fil(irow,icol) = 0.5*pi; 
                end
            end
        end

        % incrementally increase the radias from the spiral centre, until 
        % the total count of voxels excluded from the full-sized spiral 
        % exceeds a threshold 
        % threshold = expansion radius * expansion threshold
        
        % exapansion threshold: stop increment expansion of spiral radius if exceeded
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
                if centerofmass_1patt_round(1)>size(Vx_flowmap_norm_phase,2)  % remove spiral at the edge of the cortical map
                    break
                end
                if centerofmass_1patt_round(2)>size(Vx_flowmap_norm_phase,1) % remove spiral at the edge of the cortical map
                    break
                end
                if curlz(centerofmass_1patt_round(2),centerofmass_1patt_round(1),time) < 0
                spiral_filt_nega_real{ipatt,time} = sparse(-1.*temp2_01);  % if clock-wise spiral, mark it as -1 
                end
                break
            end
        end
                spiral_size_nega_real{ipatt,time} = d; % final size of each spiral at each time
      end
end


%% Similarity index (real data): find voxels within spirals that are compatible in both  phase (smoothed vs. unsmoothed) and phase vector angles (rotational)


location_matrix_x_real = [1:size(temp1_smooth_phase_real,2)];
location_matrix_y_real = [1:size(temp1_smooth_phase_real,1)];
[X_real,Y_real] = meshgrid(location_matrix_x_real,location_matrix_y_real);
cd(main_folder)
load('parcellation_template.mat')
X_real = parcellation_template(1:175,:)./parcellation_template(1:175,:).*X_real;
Y_real = parcellation_template(1:175,:)./parcellation_template(1:175,:).*Y_real;


% phase compatibility: bandpassed (smoothed) vs raw (unsmoothed)
phase_df_raw_vs_smooth_real = [];
phase_df_raw_vs_smooth_sur = [];
temp1_raw_phase_sur = temp1_raw_phase_sur.*parcellation_template(1:175,:)./parcellation_template(1:175,:);
temp1_smooth_phase_sur = temp1_smooth_phase_sur.*parcellation_template(1:175,:)./parcellation_template(1:175,:);
for t = 1:size(temp1_raw_phase_sur,3)
   temp1 = temp1_raw_phase_real(:,:,t);
   temp2 = temp1_smooth_phase_real(:,:,t);
   temp3 = temp1_raw_phase_sur(:,:,t);
   temp4 = temp1_smooth_phase_sur(:,:,t);    

   phase_df_raw_vs_smooth_real(:,:,t) = anglesubtract(temp1,temp2);
   phase_df_raw_vs_smooth_sur(:,:,t) = anglesubtract(temp3,temp4);
   
end

phase_df_raw_vs_smooth_real_abs = abs(phase_df_raw_vs_smooth_real);
temp1_raw_phase_real_filt = temp1_raw_phase_real;
temp1_raw_phase_real_filt(phase_df_raw_vs_smooth_real_abs>phase_df_threshold) = nan;

phase_df_raw_vs_smooth_sur_abs = abs(phase_df_raw_vs_smooth_sur);
temp1_raw_phase_sur_filt = temp1_raw_phase_sur;
temp1_raw_phase_sur_filt(phase_df_raw_vs_smooth_sur_abs>phase_df_threshold) = nan;

% spiral detection: phase vector (gradient) field
angle_dif_vs_distance_pos = [];
angle_dif_vs_distance_nega = [];
compatible_angle_phase_vs_distance_pos = [];
compatible_angle_phase_vs_distance_nega = [];
phaseSig = temp1_smooth_phase_real;
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

% spiral detection: anticlockwise only
curlz = [];
cav = [];
for time = 1:size(Vx_flowmap_norm_phase,3)     
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

params.minPattTime = 1; % minimum duration of spiral (10 time steps)
params.minPattSize = 3;  % minimum size of spiral (3x3 matrix)

[WCentroids,absoluteTime,instantTotalPower,pattSize,patternIdx] = pattDetection_v4(cav_filt_pos,cav_filt_pos,params,0,'CenterOfMass_Amp');
clearvars instantTotalPower patternIdx pattSize temp1_vx temp1_vy   
 
% find the centre ONLY positions in x/y coordinates
for ipatt = 1:size(WCentroids,2)
    temp1_absoluteTime = absoluteTime{ipatt};
    temp1_xy = WCentroids{ipatt};
    temp1_x = round(temp1_xy(:,1),2);
    temp1_y = round(temp1_xy(:,2),2);
    for time = 1:size(temp1_absoluteTime,2)
        spiral_filt_pos_centreONLY_real{ipatt,temp1_absoluteTime(time)} = [temp1_x(time) temp1_y(time)];
    end
end
        
% find the angle difference between center-originated vector (center of mass point to actual vector location) and phase vector
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
        temp1_vxy_angle = angle(temp1_vx + i.*temp1_vy); 
        
        % centre of mass position of the defined vector at the definded
        % time
        temp1_centerofmass_1patt = WCentroids{ipatt};
        centerofmass_1patt = temp1_centerofmass_1patt(time_withinpatt,:);
        
        % calculate the angles of spiral-centre-originated vectors toward all 
        % voxel positions across the flattened cortex
        for irow = 1:size(Vx_flowmap_norm_phase,1) 
            for icol = 1:size(Vx_flowmap_norm_phase,2) 
                y_tempalte(irow,:) = irow;
                x_tempalte(:,icol) = icol;   
            end
        end
        centerofmass_1patt_vx = round(x_tempalte - centerofmass_1patt(1));
        centerofmass_1patt_vy = round(y_tempalte - centerofmass_1patt(2));
        centerofmass_1patt_vxy_angle = angle(centerofmass_1patt_vx + i.*centerofmass_1patt_vy);
        
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
        
        angle_threshold_low = pi./2 - 0.7854; %***********%45degrees
        angle_threshold_high = pi./2 + 0.7854; %*********45degrees
        centerofmass_1patt_vxy_angle_dif_fil = centerofmass_1patt_vxy_angle_dif;
        centerofmass_1patt_vxy_angle_dif_fil(centerofmass_1patt_vxy_angle_dif_fil>angle_threshold_high) = nan;
        centerofmass_1patt_vxy_angle_dif_fil(centerofmass_1patt_vxy_angle_dif_fil<angle_threshold_low) = nan;

        centerofmass_1patt_vxy_angle_dif_fil_nonan = centerofmass_1patt_vxy_angle_dif_fil(~isnan(centerofmass_1patt_vxy_angle_dif_fil)); % voxels with compaible phase vector angles (rotational)
        distance_matrix = sqrt((X_real-centerofmass_1patt(1)).^2 + (Y_real-centerofmass_1patt(2)).^2);
        distance_matrix_nonan = distance_matrix(~isnan(centerofmass_1patt_vxy_angle_dif_fil));

        distance_matrix(isnan(distance_matrix)) = 0;
        distance_matrix_nonan_2d = zeros(size(distance_matrix));;
        distance_matrix_nonan_2d(~isnan(centerofmass_1patt_vxy_angle_dif_fil)) = distance_matrix(~isnan(centerofmass_1patt_vxy_angle_dif_fil)); % distance matrix with compatible phase vectors
        
        % find the voxels that contain both compatible phase (smoothed vs.
        % unsmoothed) and phase vector angle (to spiral)
        temp1_raw_phase_real_filt_anglealign_2d = temp1_raw_phase_real_filt(:,:,time);
        temp1_raw_phase_real_filt_anglealign_2d(isnan(centerofmass_1patt_vxy_angle_dif_fil)) = nan; 
        
        bins = [0:1:180]; % a range of possible distances (in voxels) from a sprial centre
        temp2_bin_avg = []; 
        compatible_voxel_count = [];  
        total_voxel_count = [];
        for bin_id = 1:size(bins,2)-1
            temp1_bin = centerofmass_1patt_vxy_angle_dif_fil_nonan; 
            temp1_bin(distance_matrix_nonan<bins(bin_id))=nan;
            temp1_bin(distance_matrix_nonan>bins(bin_id+1))=nan;
            temp2_bin = temp1_bin(~isnan(temp1_bin)); % voxels within a distance to sprials centres that have compatible phsae vector angles (to spiral)
            temp3_distance_bin = distance_matrix; % all voxels within a distance from spiral centres
            temp3_distance_bin(temp3_distance_bin<bins(bin_id))=nan;
            temp3_distance_bin(temp3_distance_bin>bins(bin_id+1))=nan;
            temp2_bin_avg(bin_id) = nansum(temp2_bin./temp2_bin)./nansum(temp3_distance_bin(:)./temp3_distance_bin(:)); % percentage of voxels within a set distance from sprial centres that have compatible phase vector angles (to spiral)
            temp3_phase_bin = temp1_raw_phase_real_filt_anglealign_2d; 
            temp3_phase_bin(distance_matrix_nonan_2d<bins(bin_id))=nan;
            temp3_phase_bin(distance_matrix_nonan_2d>bins(bin_id+1))=nan;   
            compatible_voxel_count(bin_id) = nansum(temp3_phase_bin(:)./temp3_phase_bin(:)); % total number of voxels within a set distance to spiral centres with compatible phases (smoothed vs. unsmoothed)
            total_voxel_count(bin_id) = nansum(temp3_distance_bin(:)./temp3_distance_bin(:)); % total number of voxels within a distance from spiral centres
        end
        angle_dif_vs_distance_pos{ipatt,time} = temp2_bin_avg; 
        % percentage of voxels with both compatible phase (smooth vs. unsmooth) and vector angle to spiral
        compatible_angle_phase_vs_distance_pos{ipatt,time} = [compatible_voxel_count(:) total_voxel_count(:)];  
      end
end


% spiral detection: clockwise spirals only
curlz = [];
cav = [];  
for time = 1:size(Vx_flowmap_norm_phase,3)      
    temp1_vx = Vx_flowmap_norm_phase(:,:,time);
    temp1_vy = Vy_flowmap_norm_phase(:,:,time);
    [x,y] = meshgrid(1:size(temp1_vx,2),1:size(temp1_vx,1));
    [curlz(:,:,time) cav(:,:,time)]  = curl(x,y,temp1_vx,temp1_vy);  % curl value of phase vector field    
end      
cav_filt_nega = curlz;
cav_filt_nega(cav_filt_nega>-1) = 0;     % detect only clockwise vorties
cav_filt_nega(isnan(cav_filt_nega)) = 0; % filter out voxels points with curl value > -1

% spiral core detection based on filtered curl values (requires function
% pattDetection_V4.m)

params.minPattTime = 1; % minimum duration of spiral (10 time steps)
params.minPattSize = 3;  % minimum size of spiral (3x3 matrix)

[WCentroids,absoluteTime,instantTotalPower,pattSize,patternIdx] = pattDetection_v4(cav_filt_nega,cav_filt_nega,params,0,'CenterOfMass_Amp');
clearvars instantTotalPower patternIdx pattSize temp1_vx temp1_vy   
 
% find the centre ONLY positions in x/y coordinates
for ipatt = 1:size(WCentroids,2)
    temp1_absoluteTime = absoluteTime{ipatt};
    temp1_xy = WCentroids{ipatt};
    temp1_x = round(temp1_xy(:,1),2);
    temp1_y = round(temp1_xy(:,2),2);
    for time = 1:size(temp1_absoluteTime,2)
        spiral_filt_nega_centreONLY_real{ipatt,temp1_absoluteTime(time)} = [temp1_x(time) temp1_y(time)];
    end
end

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
        temp1_vxy_angle = angle(temp1_vx + i.*temp1_vy); 
        
        % centre of mass position of the defined vector at the definded
        % time
        temp1_centerofmass_1patt = WCentroids{ipatt};
        centerofmass_1patt = temp1_centerofmass_1patt(time_withinpatt,:);
        
        % calculate the angles of spiral-centre-originated vectors toward all 
        % voxel positions across the flattened cortex
        for irow = 1:size(Vx_flowmap_norm_phase,1) 
            for icol = 1:size(Vx_flowmap_norm_phase,2) 
                y_tempalte(irow,:) = irow;
                x_tempalte(:,icol) = icol;   
            end
        end
        centerofmass_1patt_vx = round(x_tempalte - centerofmass_1patt(1));
        centerofmass_1patt_vy = round(y_tempalte - centerofmass_1patt(2));
        centerofmass_1patt_vxy_angle = angle(centerofmass_1patt_vx + i.*centerofmass_1patt_vy);
        
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
        
        centerofmass_1patt_vxy_angle_dif_fil_nonan = centerofmass_1patt_vxy_angle_dif_fil(~isnan(centerofmass_1patt_vxy_angle_dif_fil));
        distance_matrix = sqrt((X_real-centerofmass_1patt(1)).^2 + (Y_real-centerofmass_1patt(2)).^2);
        distance_matrix_nonan = distance_matrix(~isnan(centerofmass_1patt_vxy_angle_dif_fil));

        distance_matrix(isnan(distance_matrix)) = 0;
        distance_matrix_nonan_2d = zeros(size(distance_matrix));;
        distance_matrix_nonan_2d(~isnan(centerofmass_1patt_vxy_angle_dif_fil)) = distance_matrix(~isnan(centerofmass_1patt_vxy_angle_dif_fil));
        % find the voxels that contain both compatible phase (smoothed vs.
        % unsmoothed) and phase vector angle (rotational)
        temp1_raw_phase_real_filt_anglealign_2d = temp1_raw_phase_real_filt(:,:,time);
        temp1_raw_phase_real_filt_anglealign_2d(isnan(centerofmass_1patt_vxy_angle_dif_fil)) = nan; 
        
        
        bins = [0:1:180]; % a range of possible distances (in voxels) from a sprial centre
        temp2_bin_avg = [];
        compatible_voxel_count = [];
        total_voxel_count = [];        
        for bin_id = 1:size(bins,2)-1
            temp1_bin = centerofmass_1patt_vxy_angle_dif_fil_nonan;
            temp1_bin(distance_matrix_nonan<bins(bin_id))=nan;
            temp1_bin(distance_matrix_nonan>bins(bin_id+1))=nan;
            temp2_bin = temp1_bin(~isnan(temp1_bin));
            temp3_distance_bin = distance_matrix;
            temp3_distance_bin(temp3_distance_bin<bins(bin_id))=nan;
            temp3_distance_bin(temp3_distance_bin>bins(bin_id+1))=nan;
            temp2_bin_avg(bin_id) = nansum(temp2_bin./temp2_bin)./nansum(temp3_distance_bin(:)./temp3_distance_bin(:));
            temp3_phase_bin = temp1_raw_phase_real_filt_anglealign_2d;
            temp3_phase_bin(distance_matrix_nonan_2d<bins(bin_id))=nan;
            temp3_phase_bin(distance_matrix_nonan_2d>bins(bin_id+1))=nan;   
            compatible_voxel_count(bin_id) = nansum(temp3_phase_bin(:)./temp3_phase_bin(:));
            total_voxel_count(bin_id) = nansum(temp3_distance_bin(:)./temp3_distance_bin(:));
        end
        angle_dif_vs_distance_nega{ipatt,time} = temp2_bin_avg;
        % percentage of voxels with both compatible phase (smooth vs. unsmooth) and vector angle to spiral
        compatible_angle_phase_vs_distance_nega{ipatt,time} = [compatible_voxel_count(:) total_voxel_count(:)]; 
      end
end

% similarity index: combine the ratio of compatible voxels  within sprials from both the clockwise (nega) and anticlockwise (pos) spirals 
compatible_angle_phase_vs_distance_accu = [];
for ipatt = 1:size(compatible_angle_phase_vs_distance_nega,1);
    for t = 1:size(compatible_angle_phase_vs_distance_nega,2);
        temp1 = compatible_angle_phase_vs_distance_nega{ipatt,t};
        if nansum(temp1(:))~=0
           compatible_angle_phase_vs_distance = temp1(:,1)./temp1(:,2);           
           compatible_angle_phase_vs_distance_accu = [compatible_angle_phase_vs_distance_accu compatible_angle_phase_vs_distance(:)];    
        end
        
    end
end
for ipatt = 1:size(compatible_angle_phase_vs_distance_pos,1);
    for t = 1:size(compatible_angle_phase_vs_distance_pos,2);
        temp1 = compatible_angle_phase_vs_distance_pos{ipatt,t};
        if nansum(temp1(:))~=0
           compatible_angle_phase_vs_distance = temp1(:,1)./temp1(:,2);                       
           compatible_angle_phase_vs_distance_accu = [compatible_angle_phase_vs_distance_accu compatible_angle_phase_vs_distance(:)];    
        end
        
    end
end

compatible_angle_phase_vs_distance_accu_avg = nanmean(compatible_angle_phase_vs_distance_accu(:));


%% Similarity index (surogate data): find voxels within spirals that are compatible in both  phase (smoothed vs. unsmoothed) and phase vector angles (rotational)


location_matrix_x_sur = [1:size(temp1_smooth_phase_sur,2)];
location_matrix_y_sur = [1:size(temp1_smooth_phase_sur,1)];
[X_sur,Y_sur] = meshgrid(location_matrix_x_sur,location_matrix_y_sur);

X_sur = parcellation_template(1:175,:)./parcellation_template(1:175,:).*X_sur;
Y_sur = parcellation_template(1:175,:)./parcellation_template(1:175,:).*Y_sur;


% surrogate phase vector field

angle_dif_vs_distance_pos_sur = [];
angle_dif_vs_distance_nega_sur = [];
compatible_angle_phase_vs_distance_pos_sur = [];
compatible_angle_phase_vs_distance_nega_sur = [];
phaseSig = temp1_smooth_phase_sur;
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

% spiral detection: anticlockwise only
curlz = [];
cav = [];
for time = 1:size(Vx_flowmap_norm_phase,3)     
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

params.minPattTime = 1; % minimum duration of spiral (10 time steps)
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
        temp1_vxy_angle = angle(temp1_vx + i.*temp1_vy); 
        
        % centre of mass position of the defined vector at the definded
        % time
        temp1_centerofmass_1patt = WCentroids{ipatt};
        centerofmass_1patt = temp1_centerofmass_1patt(time_withinpatt,:);

        % calculate the angles of spiral-centre-originated vectors toward all 
        % voxel positions across the flattened cortex
        for irow = 1:size(Vx_flowmap_norm_phase,1) 
            for icol = 1:size(Vx_flowmap_norm_phase,2) 
                y_tempalte(irow,:) = irow;
                x_tempalte(:,icol) = icol;   
            end
        end
        centerofmass_1patt_vx = round(x_tempalte - centerofmass_1patt(1));
        centerofmass_1patt_vy = round(y_tempalte - centerofmass_1patt(2));
        centerofmass_1patt_vxy_angle = angle(centerofmass_1patt_vx + i.*centerofmass_1patt_vy);
        
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

        centerofmass_1patt_vxy_angle_dif_fil_nonan = centerofmass_1patt_vxy_angle_dif_fil(~isnan(centerofmass_1patt_vxy_angle_dif_fil));
        distance_matrix = sqrt((X_sur-centerofmass_1patt(1)).^2 + (Y_sur-centerofmass_1patt(2)).^2);
        distance_matrix_nonan = distance_matrix(~isnan(centerofmass_1patt_vxy_angle_dif_fil));

        distance_matrix(isnan(distance_matrix)) = 0;
        distance_matrix_nonan_2d = zeros(size(distance_matrix));;
        distance_matrix_nonan_2d(~isnan(centerofmass_1patt_vxy_angle_dif_fil)) = distance_matrix(~isnan(centerofmass_1patt_vxy_angle_dif_fil));
        temp1_raw_phase_sur_filt_anglealign_2d = temp1_raw_phase_sur_filt(:,:,time);
        temp1_raw_phase_sur_filt_anglealign_2d(isnan(centerofmass_1patt_vxy_angle_dif_fil)) = nan; 
        
        
        bins = [0:1:180];
        temp2_bin_avg = [];
        compatible_voxel_count = [];     
        total_voxel_count = [];        
        for bin_id = 1:size(bins,2)-1
            temp1_bin = centerofmass_1patt_vxy_angle_dif_fil_nonan;
            temp1_bin(distance_matrix_nonan<bins(bin_id))=nan;
            temp1_bin(distance_matrix_nonan>bins(bin_id+1))=nan;
            temp2_bin = temp1_bin(~isnan(temp1_bin));
            temp3_distance_bin = distance_matrix;
            temp3_distance_bin(temp3_distance_bin<bins(bin_id))=nan;
            temp3_distance_bin(temp3_distance_bin>bins(bin_id+1))=nan;
            temp2_bin_avg(bin_id) = nansum(temp2_bin./temp2_bin)./nansum(temp3_distance_bin(:)./temp3_distance_bin(:));
            temp3_phase_bin = temp1_raw_phase_sur_filt_anglealign_2d;
            temp3_phase_bin(distance_matrix_nonan_2d<bins(bin_id))=nan;
            temp3_phase_bin(distance_matrix_nonan_2d>bins(bin_id+1))=nan;   
            compatible_voxel_count(bin_id) = nansum(temp3_phase_bin(:)./temp3_phase_bin(:));
            total_voxel_count(bin_id) = nansum(temp3_distance_bin(:)./temp3_distance_bin(:));
        end
        angle_dif_vs_distance_pos_sur{ipatt,time} = temp2_bin_avg;
        % percentage of voxels with both compatible phase (smooth vs. unsmooth) and vector angle to spiral
        compatible_angle_phase_vs_distance_pos_sur{ipatt,time} = [compatible_voxel_count(:) total_voxel_count(:)]; 
      end
end
% spiral detection: clockwise spirals only
curlz = [];
cav = [];  
for time = 1:size(Vx_flowmap_norm_phase,3)      
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

params.minPattTime = 1; % minimum duration of spiral (10 time steps)
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
        temp1_vxy_angle = angle(temp1_vx + i.*temp1_vy); 
        
        % centre of mass position of the defined vector at the definded
        % time
        temp1_centerofmass_1patt = WCentroids{ipatt};
        centerofmass_1patt = temp1_centerofmass_1patt(time_withinpatt,:);

        % calculate the angles of spiral-centre-originated vectors toward all 
        % voxel positions across the flattened cortex
        for irow = 1:size(Vx_flowmap_norm_phase,1) 
            for icol = 1:size(Vx_flowmap_norm_phase,2) 
                y_tempalte(irow,:) = irow;
                x_tempalte(:,icol) = icol;   
            end
        end

        centerofmass_1patt_vx = round(x_tempalte - centerofmass_1patt(1));
        centerofmass_1patt_vy = round(y_tempalte - centerofmass_1patt(2));
        centerofmass_1patt_vxy_angle = angle(centerofmass_1patt_vx + i.*centerofmass_1patt_vy);
        
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
        
        centerofmass_1patt_vxy_angle_dif_fil_nonan = centerofmass_1patt_vxy_angle_dif_fil(~isnan(centerofmass_1patt_vxy_angle_dif_fil));
        distance_matrix = sqrt((X_sur-centerofmass_1patt(1)).^2 + (Y_sur-centerofmass_1patt(2)).^2);
        distance_matrix_nonan = distance_matrix(~isnan(centerofmass_1patt_vxy_angle_dif_fil));

        distance_matrix(isnan(distance_matrix)) = 0;
        distance_matrix_nonan_2d = zeros(size(distance_matrix));;
        distance_matrix_nonan_2d(~isnan(centerofmass_1patt_vxy_angle_dif_fil)) = distance_matrix(~isnan(centerofmass_1patt_vxy_angle_dif_fil));
        temp1_raw_phase_sur_filt_anglealign_2d = temp1_raw_phase_sur_filt(:,:,time);
        temp1_raw_phase_sur_filt_anglealign_2d(isnan(centerofmass_1patt_vxy_angle_dif_fil)) = nan; 
        
        bins = [0:1:180];
        temp2_bin_avg = [];
        compatible_voxel_count = [];
        total_voxel_count = [];        
        for bin_id = 1:size(bins,2)-1
            temp1_bin = centerofmass_1patt_vxy_angle_dif_fil_nonan;
            temp1_bin(distance_matrix_nonan<bins(bin_id))=nan;
            temp1_bin(distance_matrix_nonan>bins(bin_id+1))=nan;
            temp2_bin = temp1_bin(~isnan(temp1_bin));
            temp3_distance_bin = distance_matrix;
            temp3_distance_bin(temp3_distance_bin<bins(bin_id))=nan;
            temp3_distance_bin(temp3_distance_bin>bins(bin_id+1))=nan;
            temp2_bin_avg(bin_id) = nansum(temp2_bin./temp2_bin)./nansum(temp3_distance_bin(:)./temp3_distance_bin(:));
            temp3_phase_bin = temp1_raw_phase_sur_filt_anglealign_2d;
            temp3_phase_bin(distance_matrix_nonan_2d<bins(bin_id))=nan;
            temp3_phase_bin(distance_matrix_nonan_2d>bins(bin_id+1))=nan;   
            compatible_voxel_count(bin_id) = nansum(temp3_phase_bin(:)./temp3_phase_bin(:));
            total_voxel_count(bin_id) = nansum(temp3_distance_bin(:)./temp3_distance_bin(:));
        end
        angle_dif_vs_distance_nega_sur{ipatt,time} = temp2_bin_avg;
        % percentage of voxels with both compatible phase (smooth vs. unsmooth) and vector angle to spiral
        compatible_angle_phase_vs_distance_nega_sur{ipatt,time} = [compatible_voxel_count(:) total_voxel_count(:)]; 
      end
end

% similarity index: combine the ratio of compatible voxels  within sprials from both the clockwise (nega) and anticlockwise (pos) spirals 
compatible_angle_phase_vs_distance_accu_sur = [];
for ipatt = 1:size(compatible_angle_phase_vs_distance_nega_sur,1);
    for t = 1:size(compatible_angle_phase_vs_distance_nega_sur,2);
        temp1 = compatible_angle_phase_vs_distance_nega_sur{ipatt,t};
        if nansum(temp1(:))~=0
           compatible_angle_phase_vs_distance = temp1(:,1)./temp1(:,2);           
           compatible_angle_phase_vs_distance_accu_sur = [compatible_angle_phase_vs_distance_accu_sur compatible_angle_phase_vs_distance(:)];    
        end
        
    end
end
for ipatt = 1:size(compatible_angle_phase_vs_distance_pos_sur,1);
    for t = 1:size(compatible_angle_phase_vs_distance_pos_sur,2);
        temp1 = compatible_angle_phase_vs_distance_pos_sur{ipatt,t};
        if nansum(temp1(:))~=0
           compatible_angle_phase_vs_distance = temp1(:,1)./temp1(:,2);                       
           compatible_angle_phase_vs_distance_accu_sur = [compatible_angle_phase_vs_distance_accu_sur compatible_angle_phase_vs_distance(:)];    
        end        
    end
end

compatible_angle_phase_vs_distance_accu_sur_avg = nanmean(compatible_angle_phase_vs_distance_accu_sur(:));


%% select only sprials with similarity ratio (or compatible/compatibility ratio) larger than 95th percentile of the surrogate data  

% similarity index (compatible ratio) of the surrogate data across
% different distances from spiral centre (radius)
compatible_angle_phase_accu_sur = [];
total_angle_phase_accu_sur = [];
for ipatt = 1:size(compatible_angle_phase_vs_distance_nega_sur,1);
    for t = 1:size(compatible_angle_phase_vs_distance_nega_sur,2);
        temp1 = compatible_angle_phase_vs_distance_nega_sur{ipatt,t};
        if nansum(temp1(:))~=0
           compatible_angle_phase_accu_sur = [compatible_angle_phase_accu_sur temp1(:,1)];   
           total_angle_phase_accu_sur = [total_angle_phase_accu_sur temp1(:,2)];   
        end        
    end
end
for ipatt = 1:size(compatible_angle_phase_vs_distance_pos_sur,1);
    for t = 1:size(compatible_angle_phase_vs_distance_pos_sur,2);
        temp1 = compatible_angle_phase_vs_distance_pos_sur{ipatt,t};
        if nansum(temp1(:))~=0
           compatible_angle_phase_accu_sur = [compatible_angle_phase_accu_sur temp1(:,1)];    
           total_angle_phase_accu_sur = [total_angle_phase_accu_sur temp1(:,2)];   
        end
    end
end

temp1_compatible_ratio_sur = [];
temp1_compatible_ratio_sur_avg = [];
temp1_compatible_ratio_sur_std = [];
for d = 1:180 
    temp1_compatible = nansum(compatible_angle_phase_accu_sur(1:d,:),1);
    temp1_total =  nansum(total_angle_phase_accu_sur(1:d,:),1);
    temp1_compatible_ratio_sur_1d = temp1_compatible./temp1_total;
    temp1_compatible_ratio_sur_avg(d) = nanmean(temp1_compatible_ratio_sur_1d(:));
    temp1_compatible_ratio_sur_std(d) = nanstd(temp1_compatible_ratio_sur_1d(:));
    % similarity index (compatibility ratio) distributions of surrogate spirals across different radii
    temp1_compatible_ratio_sur{d} = temp1_compatible_ratio_sur_1d; 
end


% calculate similarity index of the real data and filter them using the 95th
% percentile of those from the surrogate data

% anticlockwise (pos) spirals only
temp1_compatibility_ratio_pos_real = [];
temp1_compatibility_ratio_pos_sur = [];
significant_pos = [];
total_spiral_frame_count = 0;
for ipatt = 1:size(compatible_angle_phase_vs_distance_pos,1)
    for time = 1:size(compatible_angle_phase_vs_distance_pos,2)
        temp1 = compatible_angle_phase_vs_distance_pos{ipatt,time};        
        if nansum(temp1(:)) ~=0  
           total_spiral_frame_count = total_spiral_frame_count + 1;
           % limit the similarity index (compatibility ratio) calculations within the spiral radius
           radius_expansion_limit = round(spiral_size_pos_real{ipatt,time});        
           % calculate similarity index (compatibility ratio)of the real data for each sprial at each time step
           temp1_compatibility_ratio = nansum(temp1(1:radius_expansion_limit,1))./nansum(temp1(1:radius_expansion_limit,2));
           temp1_compatibility_ratio_pos_real = [temp1_compatibility_ratio_pos_real;temp1_compatibility_ratio];
           % find the 95th percentile of similarity index (compatibility
           % ratio) from the surrogate data
           temp1_sur = temp1_compatible_ratio_sur{radius_expansion_limit};
           temp1_sur_sort = sort(temp1_sur);
           temp1_sur_sort_95perc = temp1_sur_sort(round(0.95.*size(temp1_sur_sort(:),1)));
           % select only real spiral frames with similarity index above 95th percentile of surrogate data
           if temp1_compatibility_ratio > temp1_sur_sort_95perc
              significant_pos{ipatt,time} = 1;
           end
           temp1_compatibility_ratio_pos_sur = [temp1_compatibility_ratio_pos_sur;nanmean(temp1_sur(:),1)];           
        end           
    end
end
% clockwise (nega) spirals only
temp1_compatibility_ratio_nega_real = [];
temp1_compatibility_ratio_nega_sur = [];
significant_nega = [];
for ipatt = 1:size(compatible_angle_phase_vs_distance_nega,1)
    for time = 1:size(compatible_angle_phase_vs_distance_nega,2)
        temp1 = compatible_angle_phase_vs_distance_nega{ipatt,time};
        if nansum(temp1(:)) ~=0  
           total_spiral_frame_count = total_spiral_frame_count + 1;  
           radius_expansion_limit = round(spiral_size_nega_real{ipatt,time}); % flexible expansion radius limit
           temp1_compatibility_ratio = nansum(temp1(1:radius_expansion_limit,1))./nansum(temp1(1:radius_expansion_limit,2));           
           temp1_compatibility_ratio_nega_real = [temp1_compatibility_ratio_nega_real;temp1_compatibility_ratio];
           temp1_sur = temp1_compatible_ratio_sur{radius_expansion_limit};
           temp1_sur_sort = sort(temp1_sur);
           temp1_sur_sort_95perc = temp1_sur_sort(round(0.95.*size(temp1_sur_sort(:),1)));
           if temp1_compatibility_ratio > temp1_sur_sort_95perc
              significant_nega{ipatt,time} = 1;
           end
           temp1_compatibility_ratio_nega_sur = [temp1_compatibility_ratio_nega_sur;nanmean(temp1_sur(:),1)];                   
        end                   
    end
end

% summerise the significant sprial frames for both anticlockwise and
% clockwise spirals
temp2_nega = [];
temp2_pos = [];
for ipatt_nega = 1:size(significant_nega,1)
    for time = 1:size(significant_nega,2)
        if isempty(significant_nega{ipatt_nega,time})~=1
%         temp1 = [significant_nega{ipatt_nega,time} ipatt_nega time 2];
        temp1 = [ipatt_nega time 2];
        temp2_nega = [temp2_nega;temp1 ];
        end
        
    end
end
        

for ipatt_pos = 1:size(significant_pos,1)
    for time = 1:size(significant_pos,2)
        if isempty(significant_pos{ipatt_pos,time})~=1
%         temp1 = [significant_pos{ipatt_pos,time} ipatt_pos time 1];
        temp1 = [ipatt_pos time 1];        
        temp2_pos = [temp2_pos;temp1 ];
        end
    end
end
temp2_pos_nega = [temp2_pos;temp2_nega];

temp1_compatibility_ratio_pos_nega_real = [temp1_compatibility_ratio_pos_real(:); temp1_compatibility_ratio_nega_real(:)];
temp1_compatibility_ratio_pos_nega_sur = [temp1_compatibility_ratio_pos_sur(:); temp1_compatibility_ratio_nega_sur(:)];
temp1_compatibility_ratio_pos_nega_real_avg = nanmean(temp1_compatibility_ratio_pos_nega_real(:));
temp1_compatibility_ratio_pos_nega_sur_avg = nanmean(temp1_compatibility_ratio_pos_nega_sur(:));

% select only spirals that have similarity index (compaibility ratio)
% larger than 95th percentile of surrogate data

spiral_filt_pos_real_95perc = [];
spiral_filt_nega_real_95perc = [];
spiral_size_nega_real_95perc = [];
spiral_size_pos_real_95perc = [];
spiral_size_real_95perc_accu = [];
significant_spiral_frame_count = 0;
for ispiral = 1:size(temp2_pos_nega,1)
    ipatt = temp2_pos_nega(ispiral,1);
    time = temp2_pos_nega(ispiral,2);
    pos_or_nega = temp2_pos_nega(ispiral,3);
%     temp1 = temp2_pos_nega(ispiral,1);
%     if temp1 == 1 % larger than 95 percentile of surrogate data
        significant_spiral_frame_count = significant_spiral_frame_count +1;
        if pos_or_nega == 1
        spiral_filt_pos_real_95perc{ipatt,time} = spiral_filt_pos_real{ipatt,time};
        spiral_size_real_95perc_accu = [spiral_size_real_95perc_accu; spiral_size_pos_real{ipatt,time}];
        spiral_size_pos_real_95perc{ipatt,time} = spiral_size_pos_real{ipatt,time};
        elseif pos_or_nega == 2
        spiral_filt_nega_real_95perc{ipatt,time} = spiral_filt_nega_real{ipatt,time};
        spiral_size_real_95perc_accu = [spiral_size_real_95perc_accu; spiral_size_nega_real{ipatt,time}];        
        spiral_size_nega_real_95perc{ipatt,time} = spiral_size_nega_real{ipatt,time};
        end
%     end
end

% extend spiral frames selected by including those between significant spiral frames
significant_spiral_duration_extend = [];
significant_spiral_frame_count_extend = 0;
% anticlockwise sprials only
spiral_filt_pos_real_95perc_extend = [];
spiral_filt_pos_real_centreONLY_95perc_extend = [];
for ipatt = 1:size(spiral_filt_pos_real_95perc,1)
    start_id = 0;
    time_start = [];
    time_end = [];
    for time = 1:size(spiral_filt_pos_real_95perc,2)-1
        temp1 = spiral_filt_pos_real_95perc{ipatt,time};
        temp2 = spiral_filt_pos_real_95perc{ipatt,time+1};     
        if nansum(temp1(:))~=0 && start_id == 0 
           time_start = time;
           start_id = 1;
        end
        if nansum(temp1(:)) ~= 0 && start_id == 1 && nansum(temp2(:)) == 0
           time_end = time;            
        end           
    end
    if nansum(time_start(:)) ~= 0 && nansum(time_end(:)) ~= 0
        spiral_duration_start_end = [time_start:time_end];
        for t2 = 1:size(spiral_duration_start_end(:),1)
            t3 = spiral_duration_start_end(t2);        
            spiral_filt_pos_real_95perc_extend{ipatt,t3} = spiral_filt_pos_real{ipatt,t3};
            spiral_filt_pos_real_centreONLY_95perc_extend{ipatt,t3} = spiral_filt_pos_centreONLY_real{ipatt,t3};        
        end    
        significant_spiral_frame_count_extend = significant_spiral_frame_count_extend + size(spiral_duration_start_end(:),1);
        significant_spiral_duration_extend = [significant_spiral_duration_extend size(spiral_duration_start_end(:),1)];
    end
end
% clockwise sprials only
spiral_filt_nega_real_95perc_extend = [];
spiral_filt_nega_real_centreONLY_95perc_extend = [];
for ipatt = 1:size(spiral_filt_nega_real_95perc,1)
    start_id = 0;
    time_start = [];
    time_end = [];    
    for time = 1:size(spiral_filt_nega_real_95perc,2)-1
        temp1 = spiral_filt_nega_real_95perc{ipatt,time};
        temp2 = spiral_filt_nega_real_95perc{ipatt,time+1};     
        if nansum(temp1(:))~=0 && start_id == 0 
           time_start = time;
           start_id = 1;
        end
        if nansum(temp1(:)) ~= 0 && start_id == 1 && nansum(temp2(:)) == 0
           time_end = time;            
        end           
    end
    if nansum(time_start(:)) ~= 0 && nansum(time_end(:)) ~= 0
        spiral_duration_start_end = [time_start:time_end];
        for t2 = 1:size(spiral_duration_start_end(:),1)
            t3 = spiral_duration_start_end(t2);
            spiral_filt_nega_real_95perc_extend{ipatt,t3} = spiral_filt_nega_real{ipatt,t3};
            spiral_filt_nega_real_centreONLY_95perc_extend{ipatt,t3} = spiral_filt_nega_centreONLY_real{ipatt,t3};
        end    
        significant_spiral_frame_count_extend = significant_spiral_frame_count_extend + size(spiral_duration_start_end(:),1);
        significant_spiral_duration_extend = [significant_spiral_duration_extend size(spiral_duration_start_end(:),1)];
    end
end

% radius of significant spiral frames selected for further processing
significant_spiral_radius_avg = nanmean(spiral_size_real_95perc_accu(:));
significant_spiral_duration_extend_avg = nanmean(significant_spiral_duration_extend(:));
% ratio of significant spiral frame selected for further processing
% as a percentage to all spiral frame examined
significant_spiral_frame_kept_ratio = significant_spiral_frame_count./total_spiral_frame_count;
significant_spiral_frame_kept_ratio_extend = significant_spiral_frame_count_extend./total_spiral_frame_count

if flagRest == 0 && flagTask ~= 3
    spiral_template_filt = zeros(175,251,315);
elseif flagRest == 1
    spiral_template_filt = zeros(175,251,1199);
elseif flagRest == 0 && flagTask == 3
    spiral_template_filt = zeros(175,251,405);
end
for ipatt = 1:size(spiral_filt_nega_real_95perc_extend,1)
    for t = 1:size(spiral_filt_nega_real_95perc_extend,2)
        temp1 = full(spiral_filt_nega_real_95perc_extend{ipatt,t});
        if nansum(temp1(:))~=0
           spiral_template_filt(:,:,t) = spiral_template_filt(:,:,t) + temp1;
        end
    end
end
for ipatt = 1:size(spiral_filt_pos_real_95perc_extend,1)
    for t = 1:size(spiral_filt_pos_real_95perc_extend,2)
        temp1 = full(spiral_filt_pos_real_95perc_extend{ipatt,t});
        if nansum(temp1(:))~=0
           spiral_template_filt(:,:,t) = spiral_template_filt(:,:,t) + temp1;
        end
    end
end

if flagRest == 0 
    if flagTask == 1
        if hemisphere == 1
           folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Spiral Detected'];
           cd(folder_name)
           filename = ['Spiral_detected_surfilt_language_task_orig100_LEFT_sub',num2str(subject),'.mat'];
           save([folder_name,'/',filename],'spiral_filt_nega_centreONLY_real','spiral_filt_pos_centreONLY_real','spiral_filt_nega_real','spiral_filt_pos_real','spiral_filt_nega_real_95perc_extend','spiral_filt_pos_real_95perc_extend','spiral_filt_nega_real_centreONLY_95perc_extend','spiral_filt_pos_real_centreONLY_95perc_extend','spiral_template_filt','spiral_size_real_95perc_accu','significant_spiral_radius_avg','significant_spiral_duration_extend','significant_spiral_duration_extend_avg','significant_spiral_frame_kept_ratio','significant_spiral_frame_kept_ratio_extend','expansion_threshold','temp1_compatibility_ratio_pos_nega_real','temp1_compatibility_ratio_pos_nega_sur','temp1_compatibility_ratio_pos_nega_real_avg','temp1_compatibility_ratio_pos_nega_sur_avg','spiral_filt_nega_sur','spiral_filt_pos_sur','-v7.3');
        elseif hemisphere == 2
           folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Spiral Detected'];
           cd(folder_name)
           filename = ['Spiral_detected_surfilt_language_task_orig100_RIGHT_sub',num2str(subject),'.mat'];
           save([folder_name,'/',filename],'spiral_filt_nega_centreONLY_real','spiral_filt_pos_centreONLY_real','spiral_filt_nega_real','spiral_filt_pos_real','spiral_filt_nega_real_95perc_extend','spiral_filt_pos_real_95perc_extend','spiral_filt_nega_real_centreONLY_95perc_extend','spiral_filt_pos_real_centreONLY_95perc_extend','spiral_template_filt','spiral_size_real_95perc_accu','significant_spiral_radius_avg','significant_spiral_duration_extend','significant_spiral_duration_extend_avg','significant_spiral_frame_kept_ratio','significant_spiral_frame_kept_ratio_extend','expansion_threshold','temp1_compatibility_ratio_pos_nega_real','temp1_compatibility_ratio_pos_nega_sur','temp1_compatibility_ratio_pos_nega_real_avg','temp1_compatibility_ratio_pos_nega_sur_avg','spiral_filt_nega_sur','spiral_filt_pos_sur','-v7.3');                
        end
    elseif flagTask == 2
        if hemisphere == 1
           folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Spiral Detected'];
           cd(folder_name)
           filename = ['Spiral_detected_surfilt_language_task_add100_LEFT_sub',num2str(subject),'.mat'];
           save([folder_name,'/',filename],'spiral_filt_nega_centreONLY_real','spiral_filt_pos_centreONLY_real','spiral_filt_nega_real','spiral_filt_pos_real','spiral_filt_nega_real_95perc_extend','spiral_filt_pos_real_95perc_extend','spiral_filt_nega_real_centreONLY_95perc_extend','spiral_filt_pos_real_centreONLY_95perc_extend','spiral_template_filt','spiral_size_real_95perc_accu','significant_spiral_radius_avg','significant_spiral_duration_extend','significant_spiral_duration_extend_avg','significant_spiral_frame_kept_ratio','significant_spiral_frame_kept_ratio_extend','expansion_threshold','temp1_compatibility_ratio_pos_nega_real','temp1_compatibility_ratio_pos_nega_sur','temp1_compatibility_ratio_pos_nega_real_avg','temp1_compatibility_ratio_pos_nega_sur_avg','spiral_filt_nega_sur','spiral_filt_pos_sur','-v7.3');
        elseif hemisphere == 2
           folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Spiral Detected'];
           cd(folder_name)
           filename = ['Spiral_detected_surfilt_language_task_add100_RIGHT_sub',num2str(subject),'.mat'];
           save([folder_name,'/',filename],'spiral_filt_nega_centreONLY_real','spiral_filt_pos_centreONLY_real','spiral_filt_nega_real','spiral_filt_pos_real','spiral_filt_nega_real_95perc_extend','spiral_filt_pos_real_95perc_extend','spiral_filt_nega_real_centreONLY_95perc_extend','spiral_filt_pos_real_centreONLY_95perc_extend','spiral_template_filt','spiral_size_real_95perc_accu','significant_spiral_radius_avg','significant_spiral_duration_extend','significant_spiral_duration_extend_avg','significant_spiral_frame_kept_ratio','significant_spiral_frame_kept_ratio_extend','expansion_threshold','temp1_compatibility_ratio_pos_nega_real','temp1_compatibility_ratio_pos_nega_sur','temp1_compatibility_ratio_pos_nega_real_avg','temp1_compatibility_ratio_pos_nega_sur_avg','spiral_filt_nega_sur','spiral_filt_pos_sur','-v7.3');                
        end
    elseif flagTask == 3
        if hemisphere == 1
           folder_name = [main_folder,'/Sample Data/Working Memory Task/Spiral Detected'];
           cd(folder_name)
           filename = ['Spiral_detected_surfilt_WM_task_LEFT_sub',num2str(subject),'.mat'];
           save([folder_name,'/',filename],'spiral_filt_nega_centreONLY_real','spiral_filt_pos_centreONLY_real','spiral_filt_nega_real','spiral_filt_pos_real','spiral_filt_nega_real_95perc_extend','spiral_filt_pos_real_95perc_extend','spiral_filt_nega_real_centreONLY_95perc_extend','spiral_filt_pos_real_centreONLY_95perc_extend','spiral_template_filt','spiral_size_real_95perc_accu','significant_spiral_radius_avg','significant_spiral_duration_extend','significant_spiral_duration_extend_avg','significant_spiral_frame_kept_ratio','significant_spiral_frame_kept_ratio_extend','expansion_threshold','temp1_compatibility_ratio_pos_nega_real','temp1_compatibility_ratio_pos_nega_sur','temp1_compatibility_ratio_pos_nega_real_avg','temp1_compatibility_ratio_pos_nega_sur_avg','spiral_filt_nega_sur','spiral_filt_pos_sur','-v7.3');
        elseif hemisphere == 2
           folder_name = [main_folder,'/Sample Data/Working Memory Task/Spiral Detected'];
           cd(folder_name)
           filename = ['Spiral_detected_surfilt_WM_task_RIGHT_sub',num2str(subject),'.mat'];
           save([folder_name,'/',filename],'spiral_filt_nega_centreONLY_real','spiral_filt_pos_centreONLY_real','spiral_filt_nega_real','spiral_filt_pos_real','spiral_filt_nega_real_95perc_extend','spiral_filt_pos_real_95perc_extend','spiral_filt_nega_real_centreONLY_95perc_extend','spiral_filt_pos_real_centreONLY_95perc_extend','spiral_template_filt','spiral_size_real_95perc_accu','significant_spiral_radius_avg','significant_spiral_duration_extend','significant_spiral_duration_extend_avg','significant_spiral_frame_kept_ratio','significant_spiral_frame_kept_ratio_extend','expansion_threshold','temp1_compatibility_ratio_pos_nega_real','temp1_compatibility_ratio_pos_nega_sur','temp1_compatibility_ratio_pos_nega_real_avg','temp1_compatibility_ratio_pos_nega_sur_avg','spiral_filt_nega_sur','spiral_filt_pos_sur','-v7.3');                
        end          
    end
end        
    
end
