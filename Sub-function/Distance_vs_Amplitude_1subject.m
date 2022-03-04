
function [distance_from_centre_amplitude_data_nega_accu,distance_from_centre_amplitude_data_pos_accu] = Distance_vs_Amplitude_1subject(subject)
%%  Distance from vortex centre vs fMRI amplitude: Interpolated centre
disp(['initiating process...'])    
disp(['initating processing of Distance and amplitude for single subject...'])
disp(['starting subject ',num2str(subject)])    

%% load vortex centre data
disp(['loading preprocessed fMRI data files, vortex centre coordinates and full-sized vortex distributions...'])    

filename = ['vortex_centreonly_bandpass_29.35_14.93_subject_',num2str(subject),'.mat'];
load(filename)
vortex_filt_nega_centre = vortex_filt_nega;
vortex_filt_pos_centre = vortex_filt_pos;
% load vortex data
filename = ['vortex_bandpass_29.35_14.93_subject_',num2str(subject),'.mat'];
load(filename)    
% load fMRI amplitude data
filename = ['Preprocessed_bandpass_ScaleSigma_29.35_14.93_sub',num2str(subject),'.mat'];
load(filename)    

%% calculate distance from vortex centre to each vorxels within a vortex: clockwise vortex only
disp(['processing clockwise vortices...'])    

distance_from_centre_amplitude_data_nega_accu = [];

for time = 1:size(vortex_filt_nega,2)       
    for ipatt = 1:size(vortex_filt_nega,1)
        % look for x,y coordinates of each voxel in a votex at a time point
        temp1_nega = vortex_filt_nega{ipatt,time};
        if nansum(temp1_nega(:)) == 0
          continue
        end        
        [y_nega,x_nega] = find(temp1_nega); 
        
        % look for x,y coordinates of vortex centre in a vortex at a time point
        temp1_nega_centre = vortex_filt_nega_centre{ipatt,time};
        y_nega_centre = temp1_nega_centre(2);
        x_nega_centre = temp1_nega_centre(1);
        
        % calculate the distance between vortex centre and all vorxels
        % in the vortex
        xy_nega_distance_from_centre = sqrt((x_nega_centre - x_nega).^2 + (y_nega_centre - y_nega).^2);
        
        % find the absolute amplitude value of each voxel in a votex at a time point
        xy_nega_distance_from_centre_amplitude = [];
        xy_nega_distance_from_centre_amplitude_data = [];
        for pt = 1:size(y_nega,1)
           xy_nega_distance_from_centre_amplitude(pt) = abs(sigBPass(:,y_nega(pt),x_nega(pt),time));
        end
        
        % concatenate the distance-from-centre and amplitude values of each 
        % voxel into a single variable
        xy_nega_distance_from_centre_amplitude_data = [xy_nega_distance_from_centre xy_nega_distance_from_centre_amplitude']; 
        
        % accumulatively find all the distance-from-centre and amplitude
        % values from each voxel across all vortices and time
        distance_from_centre_amplitude_data_nega_accu = [distance_from_centre_amplitude_data_nega_accu ; xy_nega_distance_from_centre_amplitude_data];
    end
end

%% calculate distance from vortex centre to each vorxels within a vortex: anticlockwise vortex only
disp(['processing anticlockwise vortices...'])    


distance_from_centre_amplitude_data_pos_accu = [];
for time = 1:size(vortex_filt_pos,2)       
    for ipatt = 1:size(vortex_filt_pos,1)
        % look for x,y coordinates of each voxel in a votex at a time point
        temp1_pos = vortex_filt_pos{ipatt,time};
        if nansum(temp1_pos(:)) == 0
          continue
        end        
        [y_pos,x_pos] = find(temp1_pos); 
        
        % look for x,y coordinates of vortex centre in a vortex at a time point
        temp1_pos_centre = vortex_filt_pos_centre{ipatt,time};
        y_pos_centre = temp1_pos_centre(2);
        x_pos_centre = temp1_pos_centre(1);
        
        % calculate the distance between vortex centre and all vorxels
        % in the vortex
        xy_pos_distance_from_centre = sqrt((x_pos_centre - x_pos).^2 + (y_pos_centre - y_pos).^2);
        
        % find the absolute amplitude value of each voxel in a votex at a time point
        xy_pos_distance_from_centre_amplitude = [];
        xy_pos_distance_from_centre_amplitude_data = [];
        for pt = 1:size(y_pos,1)
           xy_pos_distance_from_centre_amplitude(pt) = abs(sigBPass(:,y_pos(pt),x_pos(pt),time));
        end
        
        % concatenate the distance-from-centre and amplitude values of each 
        % voxel into a single variable
        xy_pos_distance_from_centre_amplitude_data = [xy_pos_distance_from_centre xy_pos_distance_from_centre_amplitude']; 
        
        % accumulatively find all the distance-from-centre and amplitude
        % values from each voxel across all vortices and time
        distance_from_centre_amplitude_data_pos_accu = [distance_from_centre_amplitude_data_pos_accu ; xy_pos_distance_from_centre_amplitude_data];
    end
end

%% save data

save_folder = [pwd,'/Sample Data/'];
save([save_folder,'Distance_vs_amplitude_bandpass_29.35_14.93_subject_',num2str(subject),'.mat'],'distance_from_centre_amplitude_data_pos_accu','distance_from_centre_amplitude_data_nega_accu')        
disp(['finishing subject ',num2str(subject)]) 
disp(['completing process...'])    

end