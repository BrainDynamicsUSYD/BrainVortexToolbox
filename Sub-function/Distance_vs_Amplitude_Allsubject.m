
function [unique_distanceidx_mm,amplitude_unique_distanceidx_avg_100subjectavg_norm] = Distance_vs_Amplitude_Allsubject(No_of_Subject)
disp(['initiating process...'])    
disp(['initating processing of distance and amplitude for all subjects...'])
%% combine the distance-from-vortex-centre and corresponding amplitude values from each subject
disp(['combining distance and amplitude for all subjects...'])    

distance_from_centre_amplitude_data_posnega_accu = [];
for subject = [1:No_of_Subject] 
    % load distance-vs-amplitude data of each subject, measured by
    % Distance_vs_Amplitude_1subject.m
    disp(['starting subject ',num2str(subject)])    
    filename = ['Distance_vs_amplitude_bandpass_29.35_14.93_subject_',num2str(subject),'.mat'];
    load(filename)

    % assign subject ID to each row (as the 3rd column) of the 
    % distance-vs-amplitude data array, 3 columns = [distnace amplitude subject]
    subject_label = subject.*ones(size(distance_from_centre_amplitude_data_pos_accu,1),1);
    distance_from_centre_amplitude_data_pos_accu = [distance_from_centre_amplitude_data_pos_accu subject_label];
    subject_label = subject.*ones(size(distance_from_centre_amplitude_data_nega_accu,1),1);
    distance_from_centre_amplitude_data_nega_accu = [distance_from_centre_amplitude_data_nega_accu subject_label];

    % combine all subjects' distance-vs-amplitude data into a variable
    distance_from_centre_amplitude_data_posnega_accu = [distance_from_centre_amplitude_data_posnega_accu; distance_from_centre_amplitude_data_pos_accu ; distance_from_centre_amplitude_data_nega_accu];
end
    
    % take 2 decimal places for the distance between vortex centre and each
    % voxel
    distance_from_centre_amplitude_data_posnega_accu_round2 = round(distance_from_centre_amplitude_data_posnega_accu(:,1),2);
    
    % find the unique values of distances across all subjects
    unique_distanceidx = unique(distance_from_centre_amplitude_data_posnega_accu_round2);
    
    % relocate the combined distance-vs-amplitude variable into 3 seperate
    % variables
    distance = distance_from_centre_amplitude_data_posnega_accu_round2;
    amplitude = distance_from_centre_amplitude_data_posnega_accu(:,2);
    subject_label = distance_from_centre_amplitude_data_posnega_accu(:,3);
    
    amplitude_unique_distanceidx_avg_1subject = [];
for subject = 1:No_of_Subject
    
    % for each subject, find the corresponding distance and amplitude
    idx = find(subject_label==subject);
    distance_1subject = distance(idx);
    amplitude_1subject = amplitude(subject_label==subject);
    amplitude_unique_distanceidx = [];
    
    % for each subject, find the amplitude values corresponding 
    % to each unique distance 
    for idx = 1:size(unique_distanceidx,1)
        idx_distance_unique_distanceidx = find(distance_1subject == unique_distanceidx(idx));
        amplitude_unique_distanceidx{idx} =  amplitude_1subject(idx_distance_unique_distanceidx);
    end
    % accmulatively combine the data from all subjects
    for idx = 1:size(unique_distanceidx,1)
        temp1_data = amplitude_unique_distanceidx{idx};
        amplitude_unique_distanceidx_avg_1subject(idx,subject) = nanmean(temp1_data(:));
    end
end

% calculate the mean, standard deviations and standard errors of ampltidue values
% corresponding to each unique distance across all subjects
for idx = 1:size(amplitude_unique_distanceidx_avg_1subject,1)
    temp1 = amplitude_unique_distanceidx_avg_1subject(idx,:);
    amplitude_unique_distanceidx_avg_100subjectavg(idx) = nanmean(temp1(:));
    amplitude_unique_distanceidx_avg_100subjectavg_std(idx) = nanstd(temp1(:));
end
amplitude_unique_distanceidx_avg_100subjectavg_stderr = amplitude_unique_distanceidx_avg_100subjectavg_std./sqrt(No_of_Subject);

%% visualization of the distance-amplitude relationship with a line plot 
disp(['visualizing dsitance-amplitude relationship...'])    

    max = nanmax(amplitude_unique_distanceidx_avg_100subjectavg(:));
    % normalization of mean amplitude corresponding to each unique distance 
    % as a percentage of the maximum value  
    amplitude_unique_distanceidx_avg_100subjectavg_norm = amplitude_unique_distanceidx_avg_100subjectavg./max;
    
    % convert voxel-based distance to mm 
    unique_distanceidx_mm = unique_distanceidx.*0.72;
    
    % line plot with standard error 
    errlow = amplitude_unique_distanceidx_avg_100subjectavg_norm - amplitude_unique_distanceidx_avg_100subjectavg_stderr./max;;
    errhigh = amplitude_unique_distanceidx_avg_100subjectavg_norm + amplitude_unique_distanceidx_avg_100subjectavg_stderr./max;;
    figure()
    plot(unique_distanceidx_mm(5:4:end),amplitude_unique_distanceidx_avg_100subjectavg_norm(5:4:end),'k')
    hold on
    plot(unique_distanceidx_mm(5:4:end),errlow(5:4:end),'color',[0.5,0.5,0.5]);
    plot(unique_distanceidx_mm(5:4:end),errhigh(5:4:end),'color',[0.5,0.5,0.5]);
    hold off
    title(['distance from vortex centre vs. amplitude'])
%     xlim([0,9])
%     ylim([0,0.8])
    xlabel(['distance from sigularity (mm)'])
    ylabel(['normalized amplitude'])

disp(['completing process...'])    

end