
function [spiral_template_timeavg_accu_stdnorm, spiral_radius_accu_avg, spiral_duration_accu_avg, spiral_count_timeavg, spiral_transverse_speed_accu_avg] = spiral_distribution_zscore_map_speed_duration_radius_count(flagRest,flagSur,flagSmooth,flagTask,hemisphere,main_folder,No_of_Subject)

if flagRest == 0
    if flagTask == 1
        foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Spiral Detected'];
        cd(foldername)
    end
end



spiral_duration_accu = [];
spiral_radius_accu = [];
spiral_template_timeavg_accu = [];
spiral_transverse_speed_accu = [];
if flagRest == 1
    spiral_count = zeros(1199,100);
elseif flagRest == 0
    spiral_count = zeros(315,100);
end

for subject = 1:No_of_Subject
    subject
    
    spiral_template = zeros(176,251,316);
    filename = ['Spiral_detected_surfilt_language_task_orig100_LEFT_sub',num2str(subject),'.mat'] ;

    if exist(filename) ~= 2
        spiral_count(:,subject) = nan;
        continue
    else
        spiral_count(:,subject) = 0;
    end
    load(filename)
    spiral_duration_accu = [spiral_duration_accu;significant_spiral_duration_extend(:)];
    spiral_radius_accu = [spiral_radius_accu;spiral_size_real_95perc_accu(:)];
    spiral_template(1:175,1:251,1:315) = abs(spiral_template_filt);
    spiral_template_timeavg = nanmean(spiral_template,3);
    spiral_template_timeavg_accu(:,:,subject) = spiral_template_timeavg;
    % transverse speed (pos nega)
    for ipatt = 1:size(spiral_filt_nega_real_centreONLY_95perc_extend,1)
        temp1_accu_1patt = [];
        for t = 1:size(spiral_filt_nega_real_centreONLY_95perc_extend,2)
            spiral_full_distribution_math_listen_matrix_avg_end = spiral_filt_nega_real_centreONLY_95perc_extend{ipatt,t};
            if nansum(spiral_full_distribution_math_listen_matrix_avg_end(:)) ~=0
               spiral_count(t,subject) =  spiral_count(t,subject) + 1;
               temp1_accu_1patt = [temp1_accu_1patt;spiral_full_distribution_math_listen_matrix_avg_end];
            end
        end
        if nansum(temp1_accu_1patt(:)./temp1_accu_1patt(:))>1
        x_distance = temp1_accu_1patt(2:end,1)-temp1_accu_1patt(1:end-1,1);
        y_distance = temp1_accu_1patt(2:end,2)-temp1_accu_1patt(1:end-1,2);
        spiral_transverse_speed_1patt = sqrt(x_distance.^2 + y_distance.^2);
        end
        spiral_transverse_speed_accu = [spiral_transverse_speed_accu;spiral_transverse_speed_1patt(:)];

    end
    for ipatt = 1:size(spiral_filt_pos_real_centreONLY_95perc_extend,1)
        temp1_accu_1patt = [];
        for t = 1:size(spiral_filt_pos_real_centreONLY_95perc_extend,2)
            spiral_full_distribution_math_listen_matrix_avg_end = spiral_filt_pos_real_centreONLY_95perc_extend{ipatt,t};
            if nansum(spiral_full_distribution_math_listen_matrix_avg_end(:)) ~=0
               spiral_count(t,subject) =  spiral_count(t,subject) + 1;                
               temp1_accu_1patt = [temp1_accu_1patt;spiral_full_distribution_math_listen_matrix_avg_end];
            end
        end
        if nansum(temp1_accu_1patt(:)./temp1_accu_1patt(:))>1
        x_distance = temp1_accu_1patt(2:end,1)-temp1_accu_1patt(1:end-1,1);
        y_distance = temp1_accu_1patt(2:end,2)-temp1_accu_1patt(1:end-1,2);
        spiral_transverse_speed_1patt = sqrt(x_distance.^2 + y_distance.^2);
        end
        spiral_transverse_speed_accu = [spiral_transverse_speed_accu;spiral_transverse_speed_1patt(:)];
    end
        
end

spiral_transverse_speed_accu_avg = nanmean(spiral_transverse_speed_accu(:));
spiral_transverse_speed_accu_std = nanstd(spiral_transverse_speed_accu(:));
spiral_transverse_speed_accu_n = size(spiral_transverse_speed_accu_avg(:),1);

spiral_template_timeavg_accu_avg = nanmean(spiral_template_timeavg_accu,3);
spiral_template_timeavg_accu_std = nanstd(spiral_template_timeavg_accu,0,3);
% spiral distribution zscore map
spiral_template_timeavg_accu_stdnorm = spiral_template_timeavg_accu_avg./spiral_template_timeavg_accu_std;

% spiral radius
spiral_radius_accu_avg = nanmean(spiral_radius_accu(:),1);
spiral_radius_accu_std = nanstd(spiral_radius_accu(:));
spiral_radius_accu_n = size(spiral_radius_accu(:),1);

% spiral duration
spiral_duration_accu_avg = nanmean(spiral_duration_accu(:),1);
spiral_duration_accu_std = nanstd(spiral_duration_accu(:));
spiral_duration_accu_n = size(spiral_duration_accu(:),1);

% spiral count
spiral_count_timeavg = nanmean(spiral_count,1);
spiral_count_subavg = nanmean(spiral_count_timeavg(:),1);
spiral_count_substd = nanstd(spiral_count_timeavg(:));
spiral_count_n = size(spiral_count_timeavg(:),1);

foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Analysis/'];
   save([foldername,'spiral_distribution_zscore_map_transversespeed_duration_radius_count.mat'],'spiral_transverse_speed_accu','spiral_transverse_speed_accu_avg','spiral_transverse_speed_accu_std','spiral_transverse_speed_accu_n',...
       'spiral_template_timeavg_accu','spiral_template_timeavg_accu_avg','spiral_template_timeavg_accu_std','spiral_template_timeavg_accu_stdnorm',...
       'spiral_radius_accu','spiral_radius_accu_avg','spiral_radius_accu_std','spiral_radius_accu_n','spiral_duration_accu','spiral_duration_accu_avg','spiral_duration_accu_std',...
       'spiral_duration_accu_n','spiral_count','spiral_count_subavg','spiral_count_substd','spiral_count_n');

cd(main_folder)


end
