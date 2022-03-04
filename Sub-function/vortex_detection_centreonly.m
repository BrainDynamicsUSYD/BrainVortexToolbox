function [vortex_filt_nega,vortex_filt_pos,subject] = vortex_detection_centreonly(subject)
%% pre-processing of fMRI data
disp(['initiating process...'])    
disp(['initating vortex centre detection...'])
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
disp(['generating phase vector field from phase map...'])

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
disp(['detecting anticlockwise brain vortex centre positions...'])    

for time = 1:size(data_allsubject,3)      
    temp1_vx = Vx_flowmap_norm_phase(:,:,time);
    temp1_vy = Vy_flowmap_norm_phase(:,:,time);
    [x,y] = meshgrid(1:size(temp1_vx,2),1:size(temp1_vx,1));
    [curlz(:,:,time) cav(:,:,time)]  = curl(x,y,temp1_vx,temp1_vy);  % curl value of phase vector field    
end      
cav_filt_pos = curlz;
cav_filt_pos(cav_filt_pos<1) = 0;     % detect only clockwise vorties
cav_filt_pos(isnan(cav_filt_pos)) = 0; % filter out voxels points with curl value > -1

% vortex core detection based on filtered curl values (requires function
% pattDetection_V4.m)

params.minPattTime = 10; % minimum duration of vortex (10 time steps)
params.minPattSize = 3;  % minimum size of vortex (3x3 matrix)

[WCentroids,absoluteTime,instantTotalPower,pattSize,patternIdx] = pattDetection_v4(cav_filt_pos,cav_filt_pos,params,0,'CenterOfMass_Amp');
clearvars instantTotalPower patternIdx pattSize temp1_vx temp1_vy   

% record anticlockwise vortex centre positions for each vortex at each time point
for ipatt = 1:size(WCentroids,2)
    temp1_absoluteTime = absoluteTime{ipatt};
    temp1_xy = WCentroids{ipatt};
    temp1_x = round(temp1_xy(:,1));
    temp1_y = round(temp1_xy(:,2));
    for time = 1:size(temp1_absoluteTime,2)
        vortex_filt_pos{ipatt,temp1_absoluteTime(time)} = [temp1_x(time) temp1_y(time)];
    end
end
    
%% Detect candidate brain vortices via curl value calcualted from phase vector field: clockwise vortex only  
% clockwise and anticlockwise vortices are detected treated seperately
disp(['detecting clockwise brain vortex centre positions...'])    


for time = 1:size(data_allsubject,3)      
    temp1_vx = Vx_flowmap_norm_phase(:,:,time);
    temp1_vy = Vy_flowmap_norm_phase(:,:,time);
    [x,y] = meshgrid(1:size(temp1_vx,2),1:size(temp1_vx,1));
    [curlz(:,:,time) cav(:,:,time)]  = curl(x,y,temp1_vx,temp1_vy);  % curl value of phase vector field    
end      
cav_filt_pos = curlz;
cav_filt_pos(cav_filt_pos>-1) = 0;     % detect only anticlockwise vorties
cav_filt_pos(isnan(cav_filt_pos)) = 0; % filter out voxels points with curl value < 1

% vortex core detection based on filtered curl values (requires function
% pattDetection_V4.m)

params.minPattTime = 10; % minimum duration of vortex (10 time steps)
params.minPattSize = 3;  % minimum size of vortex (3x3 matrix)

[WCentroids,absoluteTime,instantTotalPower,pattSize,patternIdx] = pattDetection_v4(cav_filt_pos,cav_filt_pos,params,0,'CenterOfMass_Amp');
clearvars instantTotalPower patternIdx pattSize temp1_vx temp1_vy   

% record clockwise vortex centre positions for each vortex at each time point
for ipatt = 1:size(WCentroids,2)
    temp1_absoluteTime = absoluteTime{ipatt};
    temp1_xy = WCentroids{ipatt};
    temp1_x = round(temp1_xy(:,1));
    temp1_y = round(temp1_xy(:,2));
    for time = 1:size(temp1_absoluteTime,2)
        vortex_filt_nega{ipatt,temp1_absoluteTime(time)} = [temp1_x(time) temp1_y(time)];
    end
end

%% save data

save_folder = [pwd,'/Sample Data/'];
save([save_folder,'vortex_centreonly_bandpass_29.35_14.93_subject_',num2str(subject),'.mat'],'vortex_filt_nega','vortex_filt_pos','subject')        
disp(['finishing subject ',num2str(subject)]) 

end
