function [p_value_negative_log_math_story_listen,p_value_negative_log_math_story_answer] = spiral_contrast_significance(flagRest,flagSur,hemisphere,main_folder,No_of_Subject,flagTask)

%% spiral contrast significance map: real data (math vs story listen)

% load real data spiral distribution
foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Analysis'];
cd(foldername)
filename = ['task_specific_spiral_distribution.mat'];
load(filename)

cd(main_folder)
if hemisphere == 1
    load('parcellation_template7.mat')
elseif hemisphere == 2
    load('parcellation_template22_RightBrain_subject1-100.mat')
    parcellation_template = parcellation_template22_RightBrain_100sub(:,:,1);
end

spiral_distribution_math_listen_1sub_avg = nan(175,251,No_of_Subject);
spiral_distribution_story_listen_1sub_avg = nan(175,251,No_of_Subject);
for subject = 1:No_of_Subject

% math listen mid (t=3/5)
no_of_trial = 0;
spiral_distribution_math_listen_1sub = [];
    for trial = 1:size(spiral_distribution_math_listen,1)
        temp1 = spiral_distribution_math_listen{trial,subject};
        if nansum(temp1(:))~=0
           no_of_trial = no_of_trial+ 1;
           spiral_distribution_math_listen_1sub(:,:,no_of_trial) = temp1(:,:,3);
        end       
    end
    if nansum(spiral_distribution_math_listen_1sub(:)) == 0
        continue
    end
    spiral_distribution_math_listen_1sub_avg(:,:,subject) = nanmean(spiral_distribution_math_listen_1sub,3);

spiral_distribution_story_listen_1sub = [];
no_of_trial = 0;
    for trial = 1:size(spiral_distribution_story_listen,1)
        temp1 = spiral_distribution_story_listen{trial,subject};
        if nansum(temp1(:))~=0
           no_of_trial = no_of_trial+ 1;
           spiral_distribution_story_listen_1sub(:,:,no_of_trial) = temp1(:,:,3);
        end       
    end
    if nansum(spiral_distribution_story_listen_1sub(:)) == 0
        continue
    end
    spiral_distribution_story_listen_1sub_avg(:,:,subject) = nanmean(spiral_distribution_story_listen_1sub,3);
    
end

spiral_contrast_story_math_listen_allsub_real = abs(spiral_distribution_story_listen_1sub_avg - spiral_distribution_math_listen_1sub_avg);

spiral_contrast_story_math_listen_allsub_real_avg = nanmean(spiral_contrast_story_math_listen_allsub_real,3);
spiral_contrast_story_math_listen_allsub_real_std = nanstd(spiral_contrast_story_math_listen_allsub_real,0,3);
spiral_contrast_story_math_listen_allsub_real_avg_stdnorm = spiral_contrast_story_math_listen_allsub_real_avg./spiral_contrast_story_math_listen_allsub_real_std;

%% spiral contrast significance map: real data (math vs story answer)


spiral_distribution_math_answer_1sub_avg = nan(175,251,No_of_Subject);
spiral_distribution_story_answer_1sub_avg = nan(175,251,No_of_Subject);
for subject = 1:No_of_Subject
% math answer mid (t=3/5)
no_of_trial = 0;
spiral_distribution_math_answer_1sub = [];
    for trial = 1:size(spiral_distribution_math_answer,1)
        temp1 = spiral_distribution_math_answer{trial,subject};
        if nansum(temp1(:))~=0
           no_of_trial = no_of_trial+ 1;
           spiral_distribution_math_answer_1sub(:,:,no_of_trial) = temp1(:,:,5);
        end       
    end
    if nansum(spiral_distribution_math_answer_1sub(:)) == 0
        continue
    end
    spiral_distribution_math_answer_1sub_avg(:,:,subject) = nanmean(spiral_distribution_math_answer_1sub,3);

spiral_distribution_story_answer_1sub = [];
no_of_trial = 0;
    for trial = 1:size(spiral_distribution_story_answer,1)
        temp1 = spiral_distribution_story_answer{trial,subject};
        if nansum(temp1(:))~=0
           no_of_trial = no_of_trial+ 1;
           spiral_distribution_story_answer_1sub(:,:,no_of_trial) = temp1(:,:,5);
        end       
    end
    if nansum(spiral_distribution_story_answer_1sub(:)) == 0
        continue
    end
    spiral_distribution_story_answer_1sub_avg(:,:,subject) = nanmean(spiral_distribution_story_answer_1sub,3);
end

spiral_contrast_story_math_answer_allsub_real = abs(spiral_distribution_story_answer_1sub_avg - spiral_distribution_math_answer_1sub_avg);

spiral_contrast_story_math_answer_allsub_real_avg = nanmean(spiral_contrast_story_math_answer_allsub_real,3);
spiral_contrast_story_math_answer_allsub_real_std = nanstd(spiral_contrast_story_math_answer_allsub_real,0,3);
spiral_contrast_story_math_answer_allsub_real_avg_stdnorm = spiral_contrast_story_math_answer_allsub_real_avg./spiral_contrast_story_math_answer_allsub_real_std;


%% spiral contrast significance map: surrogate data (math vs story listen)

foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Analysis'];
cd(foldername)
filename = ['task_specific_spiral_distribution_sur.mat'];
load(filename)


spiral_distribution_math_listen_1sub_avg = nan(175,251,No_of_Subject);
spiral_distribution_story_listen_1sub_avg = nan(175,251,No_of_Subject);

for subject = 1:No_of_Subject
% math listen mid (t=3/5)
no_of_trial = 0;
spiral_distribution_math_listen_1sub = [];
    for trial = 1:size(spiral_distribution_math_listen,1)
        temp1 = spiral_distribution_math_listen{trial,subject};
        if nansum(temp1(:))~=0
           no_of_trial = no_of_trial+ 1;
           spiral_distribution_math_listen_1sub(:,:,no_of_trial) = temp1(:,:,3);
           
        end       
    end
    if nansum(spiral_distribution_math_listen_1sub(:)) == 0
        continue
    end
    spiral_distribution_math_listen_1sub_avg(:,:,subject) = nanmean(spiral_distribution_math_listen_1sub,3);

no_of_trial = 0;
spiral_distribution_story_listen_1sub = [];
    for trial = 1:size(spiral_distribution_story_listen,1)
        temp1 = spiral_distribution_story_listen{trial,subject};
        if nansum(temp1(:))~=0
           no_of_trial = no_of_trial+ 1;
           spiral_distribution_story_listen_1sub(:,:,no_of_trial) = temp1(:,:,3);
           
        end       
    end
    if nansum(spiral_distribution_story_listen_1sub(:)) == 0
        continue
    end
    spiral_distribution_story_listen_1sub_avg(:,:,subject) = nanmean(spiral_distribution_story_listen_1sub,3);
    
end
spiral_contrast_story_math_listen_allsub_sur = abs(spiral_distribution_story_listen_1sub_avg - spiral_distribution_math_listen_1sub_avg);
spiral_contrast_story_math_listen_allsub_sur_avg = nanmean(spiral_contrast_story_math_listen_allsub_sur,3);

%% spiral contrast significance map: surrogate data (math vs story answer)


spiral_distribution_math_answer_1sub_avg = nan(175,251,No_of_Subject);
spiral_distribution_story_answer_1sub_avg = nan(175,251,No_of_Subject);

for subject = 1:No_of_Subject
% math answer mid (t=3/5)
no_of_trial = 0;
spiral_distribution_math_answer_1sub = [];
    for trial = 1:size(spiral_distribution_math_answer,1)
        temp1 = spiral_distribution_math_answer{trial,subject};
        if nansum(temp1(:))~=0
           no_of_trial = no_of_trial+ 1;
           spiral_distribution_math_answer_1sub(:,:,no_of_trial) = temp1(:,:,5);
           
        end       
    end
    if nansum(spiral_distribution_math_answer_1sub(:)) == 0
        continue
    end
    spiral_distribution_math_answer_1sub_avg(:,:,subject) = nanmean(spiral_distribution_math_answer_1sub,3);

no_of_trial = 0;
spiral_distribution_story_answer_1sub = [];
    for trial = 1:size(spiral_distribution_story_answer,1)
        temp1 = spiral_distribution_story_answer{trial,subject};
        if nansum(temp1(:))~=0
           no_of_trial = no_of_trial+ 1;
           spiral_distribution_story_answer_1sub(:,:,no_of_trial) = temp1(:,:,5);
           
        end       
    end
    if nansum(spiral_distribution_story_answer_1sub(:)) == 0
        continue
    end
    spiral_distribution_story_answer_1sub_avg(:,:,subject) = nanmean(spiral_distribution_story_answer_1sub,3);
    
end
spiral_contrast_story_math_answer_allsub_sur = abs(spiral_distribution_story_answer_1sub_avg - spiral_distribution_math_answer_1sub_avg);
spiral_contrast_story_math_answer_allsub_sur_avg = nanmean(spiral_contrast_story_math_answer_allsub_sur,3);



%% Calculate p-value between real and surrogate sprial contrast map


p_value = nan(175,251);
df = nan(175,251);
for irow = 1:size(spiral_contrast_story_math_listen_allsub_sur,1)
    for icol = 1:size(spiral_contrast_story_math_listen_allsub_sur,2)
        temp1_sur = spiral_contrast_story_math_listen_allsub_sur(irow,icol,:);
        if nansum(temp1_sur(:))~=0
            temp1_real = spiral_contrast_story_math_listen_allsub_real(irow,icol,:);
            temp1_sur_filt = temp1_sur(temp1_sur~=0);
            temp1_real_filt = temp1_real(temp1_real~=0);
            [h,p,ci,stats] = ttest2(temp1_sur_filt(:),temp1_real_filt(:),'Tail','left');
            p_value(irow,icol) = p;
            df(irow,icol) = stats.df;
        end
    end
end

p_value_negative_log_math_story_listen = -1.*log10(p_value);
p_value_negative_log_math_story_listen(p_value_negative_log_math_story_listen==inf) = nan;

p_value = nan(175,251);
df = nan(175,251);
for irow = 1:size(spiral_contrast_story_math_answer_allsub_sur,1)
    for icol = 1:size(spiral_contrast_story_math_answer_allsub_sur,2)
        temp1_sur = spiral_contrast_story_math_answer_allsub_sur(irow,icol,:);
        if nansum(temp1_sur(:))~=0
            temp1_real = spiral_contrast_story_math_answer_allsub_real(irow,icol,:);
            temp1_sur_filt = temp1_sur(temp1_sur~=0);
            temp1_real_filt = temp1_real(temp1_real~=0);
            [h,p,ci,stats] = ttest2(temp1_sur_filt(:),temp1_real_filt(:),'Tail','left');
            p_value(irow,icol) = p;
            df(irow,icol) = stats.df;
        end
    end
end

p_value_negative_log_math_story_answer = -1.*log10(p_value);
p_value_negative_log_math_story_answer(p_value_negative_log_math_story_answer==inf) = nan;



foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Analysis/'];
filename = ['spiral_contrast_significance.mat'];
save([foldername,filename],'p_value_negative_log_math_story_listen','p_value_negative_log_math_story_answer')

figure()
subplot(1,2,1)
pcolor(p_value_negative_log_math_story_listen)
shading interp
colorbar
colormap jet
caxis([1.3,9])
hold on
 for parcellation_ID = 1:23
    parcellation_template_1par = parcellation_template7;
    parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
    parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;
    B = bwboundaries(parcellation_template_1par,'noholes');
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2),boundary(:,1),'-','linewidth',1,'color',[1,1,1])
        
    end
 end 
 hold off
subplot(1,2,2)
pcolor(p_value_negative_log_math_story_answer)
shading interp
colorbar
colormap jet
caxis([1.3,9])
hold on
 for parcellation_ID = 1:23
    parcellation_template_1par = parcellation_template7;
    parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
    parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;
    B = bwboundaries(parcellation_template_1par,'noholes');
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2),boundary(:,1),'-','linewidth',1,'color',[1,1,1])
        
    end
 end 
 hold off



%% parcellation hierachy summary: 7 parcellations (math vs story listen)

parcellations = [1,2,3,4,5,6,7];
temp1_ind_accu = [];
for p = 1:size(parcellations,2)
    temp1_ind = p_value_negative_log_math_story_listen(parcellation_template7 ==parcellations(p));  
    temp1_ind(temp1_ind==inf) = nan;
    temp1_ind_accu{p} = temp1_ind;
    Contrast_significance_listen_ParcelAvg(p) = nansum(temp1_ind(:))./size(temp1_ind,1);
    Contrast_significance_listen_ParcelAvg_95CI(p) = nanstd(temp1_ind(:))./sqrt(size(temp1_ind,1)).*1.96;
    Contrast_significance_listen_ParcelAvg_stderr(p) = nanstd(temp1_ind(:))./sqrt(size(temp1_ind,1));
end

errlow = - Contrast_significance_listen_ParcelAvg_stderr;
errhigh =  Contrast_significance_listen_ParcelAvg_stderr;
 
 x = [1:7];
figure()
bar(x,Contrast_significance_listen_ParcelAvg)
plotName = [{'VIS'},{'SMN'},{'AUD'},{'CON'},{'DAN'},{'FPN'},{'DMN'}] ;
set(gca,'xticklabel',plotName)
hold on
er = errorbar(x,Contrast_significance_listen_ParcelAvg,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
hold off
% ylim([0,0.12])
title(['spiral contrast significance in 7 networks, Story vs math listen'])

parcellations = [1,2,3,4,5,6,7];
temp1_ind_accu = [];
for p = 1:size(parcellations,2)
    temp1_ind = p_value_negative_log_math_story_answer(parcellation_template7 ==parcellations(p));  
    temp1_ind(temp1_ind==inf) = nan;
    temp1_ind_accu{p} = temp1_ind;
    Contrast_significance_answer_ParcelAvg(p) = nansum(temp1_ind(:))./size(temp1_ind,1);
    Contrast_significance_answer_ParcelAvg_95CI(p) = nanstd(temp1_ind(:))./sqrt(size(temp1_ind,1)).*1.96;
    Contrast_significance_answer_ParcelAvg_stderr(p) = nanstd(temp1_ind(:))./sqrt(size(temp1_ind,1));
end

errlow = - Contrast_significance_answer_ParcelAvg_stderr;
errhigh =  Contrast_significance_answer_ParcelAvg_stderr;
 
x = [1:7];
figure()
bar(x,Contrast_significance_answer_ParcelAvg)
plotName = [{'VIS'},{'SMN'},{'AUD'},{'CON'},{'DAN'},{'FPN'},{'DMN'}] ;
set(gca,'xticklabel',plotName)
hold on
er = errorbar(x,Contrast_significance_answer_ParcelAvg,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
hold off
% ylim([0,0.12])
title(['spiral contrast significance in 7 networks, Story vs math answer'])

foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Analysis/'];
filename = ['spiral_contrast_significance_7networks.mat'];
save([foldername,filename],'plotName','Contrast_significance_listen_ParcelAvg','Contrast_significance_listen_ParcelAvg_stderr','Contrast_significance_answer_ParcelAvg','Contrast_significance_answer_ParcelAvg_stderr')



end