function [count_repulsion_percent_avg,count_partialannihilation_percent_avg,count_fullannihilation_percent_avg] = spiral_interaction_statistics_revision(flagRest,flagSur,hemisphere,main_folder,No_of_Subject,flagTask)

%% spiral interaction tracking

% clear all
main_folder = pwd;
% no_subject = 180;
% disp('extract right brain data...')
% for interaction_type = 1:3
    repulsion_count = 0;
    full_annihilation_distance = 5;
    partial_annihilation_distance = 5;
    repulsion_approach_distance = 10; % specific 2 filters, approach and repulse distances
    repulsion_repulse_distance = 15; % 1st approach, then repulse   
    interaction_template_repulsion_allsubject = [];
    interaction_template_partialannihilation_allsubject = [];
    interaction_template_fullannihilation_allsubject = [];
for interaction_type = 1:3 % full annihilation =1, partial annihilation = 2, repulsion = 3

    for subject = 1:No_of_Subject
        disp(['interaction type ',num2str(interaction_type),', subject ',num2str(subject)])
        interaction_template = zeros(176,251);
        full_annihilation_accu = [];
        partial_annihilation_accu = [];
        repulsion_accu = [];
        foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Spiral Detected'];
        cd(foldername)
        filename = ['Spiral_detected_surfilt_language_task_orig100_LEFT_sub',num2str(subject),'.mat'];
        load(filename)
        
        spiral_centre_distribution_nega_1pattnotend = [];
        for ipatt = 1:size(spiral_filt_nega_real_centreONLY_95perc_extend,1)
            for time = 1:size(spiral_filt_nega_real_centreONLY_95perc_extend,2)-4
                temp1 = spiral_filt_nega_real_centreONLY_95perc_extend{ipatt,time};
                temp2 = spiral_filt_nega_real_centreONLY_95perc_extend{ipatt,time+1};
                % excluding 4s before end, so the partial annihilation only works if the surviving spiral survive 5s more
                temp4 = spiral_filt_nega_real_centreONLY_95perc_extend{ipatt,time+4}; % excluding 4s before end

                if nansum(temp1(:)) ==0 && nansum(temp2(:)) == 0
                    continue
                elseif nansum(temp1(:)) ~=0 && nansum(temp2(:)) ~= 0  && nansum(temp4(:)) ~= 0 % excluding 4s before end, so the partial annihilation only works if the surviving spiral survive 5s more 
                    temp3 = [temp1 time ipatt];
                    spiral_centre_distribution_nega_1pattnotend = [spiral_centre_distribution_nega_1pattnotend; temp3]; % excluding the end point
                elseif nansum(temp1(:)) ~=0 && nansum(temp2(:)) == 0
                    spiral_centre_distribution_nega_1pattend(ipatt,:) = [temp1 time];
                    break
                end
            end
        end
        spiral_centre_distribution_pos_1pattnotend = [];
        for ipatt = 1:size(spiral_filt_pos_real_centreONLY_95perc_extend,1)
            for time = 1:size(spiral_filt_pos_real_centreONLY_95perc_extend,2)-4
                temp1 = spiral_filt_pos_real_centreONLY_95perc_extend{ipatt,time};
                temp2 = spiral_filt_pos_real_centreONLY_95perc_extend{ipatt,time+1};
                % excluding 4s before end, so the partial annihilation only works if the surviving spiral survive 5s more
                temp4 = spiral_filt_pos_real_centreONLY_95perc_extend{ipatt,time+4}; % excluding 4s before end

                if nansum(temp1(:)) ==0 && nansum(temp2(:)) == 0
                    continue
                elseif nansum(temp1(:)) ~=0 && nansum(temp2(:)) ~= 0 && nansum(temp4(:)) ~= 0 % excluding 4s before end, so the partial annihilation only works if the surviving spiral survive 5s more  
                    temp3 = [temp1 time ipatt];
                    spiral_centre_distribution_pos_1pattnotend = [spiral_centre_distribution_pos_1pattnotend; temp3]; % excluding the end point   
                elseif nansum(temp1(:)) ~=0 && nansum(temp2(:)) == 0
                    spiral_centre_distribution_pos_1pattend(ipatt,:) = [temp1 time];
                    break
                end
            end
        end

        idx_notzero = find(spiral_centre_distribution_nega_1pattend(:,1) ~=0);
        spiral_centre_distribution_nega_1pattend = spiral_centre_distribution_nega_1pattend(idx_notzero,:);
        idx_notzero = find(spiral_centre_distribution_nega_1pattnotend(:,1) ~=0);
        spiral_centre_distribution_nega_1pattnotend = spiral_centre_distribution_nega_1pattnotend(idx_notzero,:);

        idx_notzero = find(spiral_centre_distribution_pos_1pattend(:,1) ~=0);
        spiral_centre_distribution_pos_1pattend = spiral_centre_distribution_pos_1pattend(idx_notzero,:);
        idx_notzero = find(spiral_centre_distribution_pos_1pattnotend(:,1) ~=0);
        spiral_centre_distribution_pos_1pattnotend = spiral_centre_distribution_pos_1pattnotend(idx_notzero,:);




        if interaction_type == 1 % full annihilation

            for time = 1:315
                temp1 = spiral_centre_distribution_nega_1pattend(:,3);
                idx = find(temp1 == time);
                spiral_centre_distribution_nega_1pattend_1time = spiral_centre_distribution_nega_1pattend(idx,:);
                temp1 = spiral_centre_distribution_pos_1pattend(:,3);
                idx = find(temp1 == time);
                spiral_centre_distribution_pos_1pattend_1time = spiral_centre_distribution_pos_1pattend(idx,:);

                 for ipatt_nega = 1:size(spiral_centre_distribution_nega_1pattend_1time,1)
                     for ipatt_pos = 1:size(spiral_centre_distribution_pos_1pattend_1time,1)
                         pos_x = spiral_centre_distribution_pos_1pattend_1time(ipatt_pos,1);
                         pos_y = spiral_centre_distribution_pos_1pattend_1time(ipatt_pos,2);
                         time = spiral_centre_distribution_pos_1pattend_1time(ipatt_pos,3);
                         nega_x = spiral_centre_distribution_nega_1pattend_1time(ipatt_nega,1);
                         nega_y = spiral_centre_distribution_nega_1pattend_1time(ipatt_nega,2);
                         distance(ipatt_nega,ipatt_pos) = sqrt((pos_x-nega_x).^2 + (pos_y-nega_y).^2); 
                         full_annhilation = [];
                         if distance(ipatt_nega,ipatt_pos) < full_annihilation_distance
                             full_annhilation = [pos_x pos_y nega_x nega_y time];
                         end
                         if nansum(full_annhilation(:)) == 0
                             continue
                         end
                         full_annihilation_accu = [full_annihilation_accu ; full_annhilation];
                     end
                 end
            end
        elseif interaction_type == 2 % partial annihilation

            for time = 1:315
                temp1 = spiral_centre_distribution_nega_1pattend(:,3);
                idx = find(temp1 == time);
                spiral_centre_distribution_nega_1pattend_1time = spiral_centre_distribution_nega_1pattend(idx,:);
                temp1 = spiral_centre_distribution_nega_1pattnotend(:,3);
                idx = find(temp1 == time);
                spiral_centre_distribution_nega_1pattnotend_1time = spiral_centre_distribution_nega_1pattnotend(idx,:);

                temp1 = spiral_centre_distribution_pos_1pattend(:,3);
                idx = find(temp1 == time);
                spiral_centre_distribution_pos_1pattend_1time = spiral_centre_distribution_pos_1pattend(idx,:);
                temp1 = spiral_centre_distribution_pos_1pattnotend(:,3);
                idx = find(temp1 == time);
                spiral_centre_distribution_pos_1pattnotend_1time = spiral_centre_distribution_pos_1pattnotend(idx,:);
                 % 1st find pairs between nega_end and pos_notend
                 distance = [];
                 for ipatt_nega = 1:size(spiral_centre_distribution_nega_1pattend_1time,1)
                     for ipatt_pos = 1:size(spiral_centre_distribution_pos_1pattnotend_1time,1)
                         pos_x = spiral_centre_distribution_pos_1pattnotend_1time(ipatt_pos,1);
                         pos_y = spiral_centre_distribution_pos_1pattnotend_1time(ipatt_pos,2);
                         time = spiral_centre_distribution_pos_1pattnotend_1time(ipatt_pos,3);
                         nega_x = spiral_centre_distribution_nega_1pattend_1time(ipatt_nega,1);
                         nega_y = spiral_centre_distribution_nega_1pattend_1time(ipatt_nega,2);
                         distance(ipatt_nega,ipatt_pos) = sqrt((pos_x-nega_x).^2 + (pos_y-nega_y).^2); 
                         partial_annhilation = [];
                         if distance(ipatt_nega,ipatt_pos) < partial_annihilation_distance % 5 distance as limit
                             partial_annhilation = [nega_x nega_y pos_x pos_y  time];
                         end
                         if nansum(partial_annhilation(:)) == 0
                             continue
                         end
                         partial_annihilation_accu = [partial_annihilation_accu ; partial_annhilation];
                     end
                 end
                 % 2nd find pairs between pos_end and nega_notend
                 distance = [];
                 for ipatt_pos = 1:size(spiral_centre_distribution_pos_1pattend_1time,1)    
                     for ipatt_nega = 1:size(spiral_centre_distribution_nega_1pattnotend_1time,1)
                         pos_x = spiral_centre_distribution_pos_1pattend_1time(ipatt_pos,1);
                         pos_y = spiral_centre_distribution_pos_1pattend_1time(ipatt_pos,2);
                         time = spiral_centre_distribution_pos_1pattend_1time(ipatt_pos,3);
                         nega_x = spiral_centre_distribution_nega_1pattnotend_1time(ipatt_nega,1);
                         nega_y = spiral_centre_distribution_nega_1pattnotend_1time(ipatt_nega,2);
                         distance(ipatt_nega,ipatt_pos) = sqrt((pos_x-nega_x).^2 + (pos_y-nega_y).^2); 
                         partial_annhilation = [];
                         if distance(ipatt_nega,ipatt_pos) < partial_annihilation_distance % 5 distance as limit
                             partial_annhilation = [pos_x pos_y nega_x nega_y time];
                         end
                         if nansum(partial_annhilation(:)) == 0
                             continue
                         end
                         partial_annihilation_accu = [partial_annihilation_accu ; partial_annhilation];
                     end
                 end

            end


        elseif interaction_type == 3 % repulsion

            for time = 1:315
                temp1 = spiral_centre_distribution_nega_1pattnotend(:,3);
                idx = find(temp1 == time);
                spiral_centre_distribution_nega_1pattnotend_1time = spiral_centre_distribution_nega_1pattnotend(idx,:);
                temp1 = spiral_centre_distribution_pos_1pattnotend(:,3);
                idx = find(temp1 == time);
                spiral_centre_distribution_pos_1pattnotend_1time = spiral_centre_distribution_pos_1pattnotend(idx,:);
                 % 1st find pairs between nega_end and pos_notend
                 distance = [];
                 for ipatt_nega = 1:size(spiral_centre_distribution_nega_1pattnotend_1time,1)
                     for ipatt_pos = 1:size(spiral_centre_distribution_pos_1pattnotend_1time,1)
                         pos_x = spiral_centre_distribution_pos_1pattnotend_1time(ipatt_pos,1);
                         pos_y = spiral_centre_distribution_pos_1pattnotend_1time(ipatt_pos,2);
                         time = spiral_centre_distribution_pos_1pattnotend_1time(ipatt_pos,3); % only need 1 time since pos and nega share the same time 
                         ipatt_idx_pos = spiral_centre_distribution_pos_1pattnotend_1time(ipatt_pos,4);                     
                         nega_x = spiral_centre_distribution_nega_1pattnotend_1time(ipatt_nega,1);
                         nega_y = spiral_centre_distribution_nega_1pattnotend_1time(ipatt_nega,2);
                         ipatt_idx_nega = spiral_centre_distribution_nega_1pattnotend_1time(ipatt_nega,4);
                         distance(ipatt_nega,ipatt_pos) = sqrt((pos_x-nega_x).^2 + (pos_y-nega_y).^2); 
                         repulsion = [];
                         if distance(ipatt_nega,ipatt_pos) < repulsion_approach_distance % 5 distance as limit
                             repulsion = [nega_x nega_y pos_x pos_y time ipatt_idx_nega ipatt_idx_pos distance(ipatt_nega,ipatt_pos)];
                         end
                         if nansum(repulsion(:)) == 0
                             continue
                         end
                         repulsion_accu = [repulsion_accu ; repulsion];
                     end
                 end

            end


        end

    % refine and consolidate results
        if interaction_type == 1
            for ipatt = 1:size(full_annihilation_accu,1)
            pos_x = round(full_annihilation_accu(ipatt,1));
            pos_y = round(full_annihilation_accu(ipatt,2));
            nega_x = round(full_annihilation_accu(ipatt,3));
            nega_y = round(full_annihilation_accu(ipatt,4));   
            interaction_template(pos_y,pos_x) = interaction_template(pos_y,pos_x) +1;
            interaction_template(nega_y,nega_x) = interaction_template(nega_y,nega_x) +1;
            interaction_template_fullannihilation = interaction_template;
            end
            interaction_template_fullannihilation_allsubject{subject} = interaction_template_fullannihilation; 
        elseif interaction_type == 2
            for ipatt = 1:size(partial_annihilation_accu,1)
            pos_x = round(partial_annihilation_accu(ipatt,1));
            pos_y = round(partial_annihilation_accu(ipatt,2));
            nega_x = round(partial_annihilation_accu(ipatt,3));
            nega_y = round(partial_annihilation_accu(ipatt,4));     
            interaction_template(pos_y,pos_x) = interaction_template(pos_y,pos_x) +1;
            interaction_template(nega_y,nega_x) = interaction_template(nega_y,nega_x) +1;
            interaction_template_partialannihilation = interaction_template;
            end
            interaction_template_partialannihilation_allsubject{subject} = interaction_template_partialannihilation;         
        elseif interaction_type == 3
            repulsion_accu_unique_repulsionconfirmed = [];
            patt_idx_past = [];
            repulsion_accu_unique_accu = [];
            for ipatt = 1:size(repulsion_accu,1)
            pos_x = round(repulsion_accu(ipatt,1));
            pos_y = round(repulsion_accu(ipatt,2));
            nega_x = round(repulsion_accu(ipatt,3));
            nega_y = round(repulsion_accu(ipatt,4));  
            time = repulsion_accu(ipatt,5);      
            patt_idx = repulsion_accu(ipatt,6:7);
            distance = repulsion_accu(ipatt,8);
            % remove duplicated interactions between the same spiral pair,
            if nansum(patt_idx) == nansum(patt_idx_past); 
                continue
            end
            patt_idx_past = patt_idx;
            repulsion_accu_unique = [nega_x nega_y pos_x pos_y time patt_idx distance ]; % find unique approach events, initial point of contact (distance <5)
            repulsion_accu_unique_accu = [repulsion_accu_unique_accu;repulsion_accu_unique];
            end
    %         find true repulsion event where distance first reduced to <5 then
    %         increased to >10 in 5s
            for ipatt = 1:size(repulsion_accu_unique_accu,1)
                patt_idx_nega = repulsion_accu_unique_accu(ipatt,6);
                patt_idx_pos = repulsion_accu_unique_accu(ipatt,7);
                time = repulsion_accu_unique_accu(ipatt,5);
                distance_raw = repulsion_accu_unique_accu(ipatt,8);
                centre_location_nega = [];
                centre_location_pos = [];
                for t = 1:4
                    temp1_nega = spiral_filt_nega_real_centreONLY_95perc_extend{patt_idx_nega,time+t};
                    temp1_pos = spiral_filt_pos_real_centreONLY_95perc_extend{patt_idx_pos,time+t};
                    %******************************* ?????
                    if nansum(temp1_nega(:)) == 0 | nansum(temp1_pos(:)) == 0 
                        continue
                    end
                    centre_location_nega(t,:) = spiral_filt_nega_real_centreONLY_95perc_extend{patt_idx_nega,time+t};
                    centre_location_pos(t,:) = spiral_filt_pos_real_centreONLY_95perc_extend{patt_idx_pos,time+t};
                end
                distance = sqrt((centre_location_nega(1)-centre_location_pos(1)).^2 + (centre_location_nega(2)-centre_location_pos(2)).^2);
                if nanmax(distance(:)) > repulsion_repulse_distance % distance filter for true repulsion 
                    repulsion_accu_unique_repulsionconfirmed = [repulsion_accu_unique_repulsionconfirmed; repulsion_accu_unique_accu(ipatt,:)] ; 
                end
            end
            if nansum(repulsion_accu_unique_repulsionconfirmed) == 0 % if no repulsion event, skip next
                continue
            elseif nansum(repulsion_accu_unique_repulsionconfirmed) ~= 0
                repulsion_count = repulsion_count + 1;
            end
            for ipatt = 1:size(repulsion_accu_unique_repulsionconfirmed,1)
            pos_x = repulsion_accu_unique_repulsionconfirmed(ipatt,1);
            pos_y = repulsion_accu_unique_repulsionconfirmed(ipatt,2);
            nega_x = repulsion_accu_unique_repulsionconfirmed(ipatt,3);
            nega_y = repulsion_accu_unique_repulsionconfirmed(ipatt,4);     
            interaction_template(pos_y,pos_x) = interaction_template(pos_y,pos_x) +1;
            interaction_template(nega_y,nega_x) = interaction_template(nega_y,nega_x) +1;
            interaction_template_repulsion = interaction_template;
            end
            interaction_template_repulsion_allsubject{subject} = interaction_template_repulsion;

        end
    end

end
%%
    count_partialannihilation = [];
    count_repulsion = [];
    count_fullannihilation = [];
for interaction_type = 1:3
    for subject = 1:No_of_Subject
        if interaction_type == 3
            interaction_template_repulsion = interaction_template_repulsion_allsubject{subject};
            count_repulsion(subject) = nansum(interaction_template_repulsion(:));
        elseif interaction_type == 2
            interaction_template_partialannihilation = interaction_template_partialannihilation_allsubject{subject};
            count_partialannihilation(subject) = nansum(interaction_template_partialannihilation(:));
        elseif   interaction_type == 1 
            interaction_template_fullannihilation = interaction_template_fullannihilation_allsubject{subject};
            count_fullannihilation(subject) = nansum(interaction_template_fullannihilation(:));
        end   

    end
end
count_total_interaction = count_repulsion + count_partialannihilation + count_fullannihilation;
count_repulsion_percent = count_repulsion./count_total_interaction;
count_repulsion_percent_avg = nanmean(count_repulsion_percent);
count_repulsion_percent_std = nanstd(count_repulsion_percent);

count_partialannihilation_percent = count_partialannihilation./count_total_interaction;
count_partialannihilation_percent_avg = nanmean(count_partialannihilation_percent(:));
count_partialannihilation_percent_std = nanstd(count_partialannihilation_percent(:));

count_fullannihilation_percent = count_fullannihilation./count_total_interaction;
count_fullannihilation_percent_avg = nanmean(count_fullannihilation_percent(:));
count_fullannihilation_percent_std = nanstd(count_fullannihilation_percent(:));

    
foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Analysis/'];
filename = ['spiral_interaction_statistics.mat'];
save([foldername,filename],'repulsion_approach_distance','repulsion_repulse_distance'...
    ,'full_annihilation_distance','partial_annihilation_distance',...
    'count_repulsion_percent_avg','count_repulsion_percent_std',...
    'count_partialannihilation_percent_avg','count_partialannihilation_percent_std'...
    ,'count_fullannihilation_percent_avg','count_fullannihilation_percent_std'...
    ,'count_fullannihilation_percent','count_partialannihilation_percent','count_repulsion_percent');
cd(main_folder)


count_interaction_typles_avg = [count_fullannihilation_percent_avg count_partialannihilation_percent_avg count_repulsion_percent_avg];
count_interaction_typles_errh = [count_fullannihilation_percent_std count_partialannihilation_percent_std count_repulsion_percent_avg+count_repulsion_percent_std];
figure()
bar(count_interaction_typles_avg,'g')
somenames={'FA'; 'PA'; 'RE'; }
set(gca,'xticklabel',somenames)
hold on
er = errorbar(1:3,count_interaction_typles_avg,count_interaction_typles_errh,count_interaction_typles_errh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
scatter(ones(size(count_fullannihilation_percent(:))).*(1+(rand(size(count_fullannihilation_percent(:)))-0.5)/3),count_fullannihilation_percent(:),'r','filled')
scatter(ones(size(count_partialannihilation_percent(:))).*(1+(rand(size(count_partialannihilation_percent(:)))+2.5)/3),count_partialannihilation_percent(:),'r','filled')
scatter(ones(size(count_repulsion_percent(:))).*(1+(rand(size(count_repulsion_percent(:)))+5.5)/3),count_repulsion_percent(:),'r','filled')
hold off
title(['vortex interaction types'])
% ylim([0,0.6])


end