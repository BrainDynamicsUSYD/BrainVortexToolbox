function [classification_accuracy_avg,classification_accuracy_stderr] = PhaseVectorAngle_local_classifier(flagSur,hemisphere,main_folder,No_of_Subject,flagTask)
disp(['initiating...'])  


foldername = [main_folder,'/Sample Data/Language Task Original 100 sub/Analysis'];
cd(foldername)
filename = ['task_specific_phase_vector_field.mat'];
load(filename) 

vx_math_listen = Vx_math_listen_accu(:,:,6:10,:);
vy_math_listen = Vy_math_listen_accu(:,:,6:10,:);

vx_math_answer = Vx_math_answer_accu(:,:,1:5,:);
vy_math_answer = Vy_math_answer_accu(:,:,1:5,:);




%% classify task conditions (math listening vs. answering) based on single-trial moment-by-moment phase vector field angles 
disp(['classify task conditions (math listening vs. answering) based on local (3x3) phase vector angles...'])    



% isolate 3x3 matrices across the Area of Interest (x:130-220, y:10-170)


classification_accuracy_avg = nan(176,251);

% for each target voxel within the area of interest (ROI,x:130-220, y:10-170), 
% find its surrounding voxels to form a 3x3 array centred at the target voxel 

for irow = 10:170
    if mod(irow,10) == 0    
    disp(['this process may take a few minutes, progress: ',num2str(round(((irow-10)./160).*100)),'%'])    
    end
    for icol = 130:220
        % find the 3x3 array centred at the target voxel
        x_array = [icol-1:icol+1];
        y_array = [irow-1:irow+1];
        % extract the angles of phase vectors within that 3x3 array
        vxy_math_listen_angle_ROI = angle(vx_math_listen(y_array,x_array,:,:) + i.*vy_math_listen(y_array,x_array,:,:));    
        vxy_math_answer_angle_ROI = angle(vx_math_answer(y_array,x_array,:,:) + i.*vy_math_answer(y_array,x_array,:,:));    

        % for each voxel, iteratively classify task conditions based on the 
        % angles of its surrounding phase vectors (3x3 array) 
        classification_accuracy = [];
        % find no. of trials for both math and story tasks
        trials_math_listen = size(vxy_math_listen_angle_ROI,4);
        trials_math_answer = size(vxy_math_answer_angle_ROI,4);

        for iteration = 1:100
        % 80% random trials for training
        train_trial_math_listen = randperm(trials_math_listen,round(0.8*trials_math_listen));
        train_trial_math_answer = randperm(trials_math_answer,round(0.8*trials_math_answer));
        % remaining 20% random trials for testing
        test_trial_math_listen = setdiff([1:trials_math_listen],train_trial_math_listen);
        test_trial_math_answer = setdiff([1:trials_math_answer],train_trial_math_answer);
        
        % calculate the mean phase vector angles from training trials 
        % for each task condition, as training template
        % each trial lasts 5 time steps
        
        % training template for math listening tasks
        temp1_math_listen_train_template_vx = nanmean(nanmean(vx_math_listen(y_array,x_array,:,:),4),3);
        temp1_math_listen_train_template_vy = nanmean(nanmean(vy_math_listen(y_array,x_array,:,:),4),3);
        temp1_math_listen_train_template_vxy_angle = angle(temp1_math_listen_train_template_vx + i.*temp1_math_listen_train_template_vy);
        % training template for math listening tasks                        
        temp1_math_answer_train_template_vx = nanmean(nanmean(vx_math_answer(y_array,x_array,:,:),4),3);
        temp1_math_answer_train_template_vy = nanmean(nanmean(vy_math_answer(y_array,x_array,:,:),4),3);
        temp1_math_answer_train_template_vxy_angle = angle(temp1_math_answer_train_template_vx + i.*temp1_math_answer_train_template_vy); 

        % extract the mean angle of phase vectors of each test trial across 5
        % time steps selected from each session
        temp1_math_listen_test_angle = permute(nanmean(vxy_math_listen_angle_ROI(:,:,:,test_trial_math_listen),3),[1,2,4,3]);
        temp1_math_answer_test_angle = permute(nanmean(vxy_math_answer_angle_ROI(:,:,:,test_trial_math_answer),3),[1,2,4,3]);

        % calculate correlation coefficient with 2 trained templates
        temp1_train_math_listen = temp1_math_listen_train_template_vxy_angle;
        temp1_train_math_answer = temp1_math_answer_train_template_vxy_angle;

        % classification of task conditions with only math LISTENING trials 
        % for each math test trial, calculate its correlation coefficient with  
        % 2 training templates (math listening and answering)
        min_trial_count = min(size(test_trial_math_listen,2),size(test_trial_math_answer,2));
        for trial = 1:min_trial_count
            % remove voxels in both test trials and training templates if 0 
            % and nan values are found in either 
            % this step is to ensure the same data size with no nan/0 in 
            % all variables when calculating correlation coefficients
            temp1_test_math_listen = temp1_math_listen_test_angle(:,:,trial);
            temp1_train_math_listen_filt = temp1_train_math_listen(~isnan(temp1_test_math_listen));
            temp1_train_math_answer_filt = temp1_train_math_answer(~isnan(temp1_test_math_listen));
            temp1_test_math_listen = temp1_test_math_listen(~isnan(temp1_test_math_listen));
            temp1_train_math_listen_filt_01 = temp1_train_math_listen_filt./temp1_train_math_listen_filt;
            temp1_train_math_answer_filt_01 = temp1_train_math_answer_filt./temp1_train_math_answer_filt;
            temp1_test_math_listen_01 = temp1_test_math_listen./temp1_test_math_listen;
            filt_notnan_01 =  temp1_train_math_listen_filt_01(:).*temp1_train_math_answer_filt_01(:).*temp1_test_math_listen_01(:);
            temp1_train_math_answer_filt_filt = temp1_train_math_answer_filt(~isnan(filt_notnan_01));
            temp1_train_math_listen_filt_filt = temp1_train_math_listen_filt(~isnan(filt_notnan_01));
            temp1_test_math_listen_filt = temp1_test_math_listen(~isnan(filt_notnan_01));
            
            if size(temp1_test_math_listen_filt(:)) == 1 % remove test trials with only 1 data points
                continue
            end
            % calcualte angle differences between phase vectors of
            % math listening test trials and math listening
            % training template
            error_listen_listen = temp1_test_math_listen_filt(:)-temp1_train_math_listen_filt_filt(:);
            % calcualte angle differences between phase vectors of
            % math listening test trials and math answering
            % training template     
            error_listen_answer = temp1_test_math_listen_filt(:)-temp1_train_math_answer_filt_filt(:);
            % restrict angle difference between -pi and pi
            for it = 1:size(error_listen_listen(:),1)
                if error_listen_listen(it)>pi
                    error_listen_listen(it) = error_listen_listen(it) - 2.*pi;
                elseif error_listen_listen(it) < -pi
                    error_listen_listen(it) = error_listen_listen(it) + 2.*pi;
                end     
            end
            for it = 1:size(error_listen_answer(:),1)
                if error_listen_answer(it)>pi
                    error_listen_answer(it) = error_listen_answer(it) - 2.*pi;
                elseif error_listen_answer(it) < -pi
                    error_listen_answer(it) = error_listen_answer(it) + 2.*pi;
                end     
            end            
            % calcualte mean-squared erros (MSE) of angle differences 
            % between phase vectors of math listening test trials
            % and math listening training template  
            MSE_listen_listen(trial) = nansum(error_listen_listen.^2)./nansum(error_listen_listen(:)./error_listen_listen(:));
            % calcualte mean-squared erros (MSE) of angle differences 
            % between phase vectors of math listening test trials
            % and math answering training template              
            MSE_listen_answer(trial) = nansum(error_listen_answer.^2)./nansum(error_listen_answer(:)./error_listen_answer(:));

        end
            % compare MSE between listen-listen and listen-answer pairs
            % (i.e., listen-listen = listening test trial vs. listening training template)
            % the one with lower MSE is selected as the predicted task condition  
            MSE_math_listen = [MSE_listen_listen' MSE_listen_answer'];
            max_count_math_listen = nan(1,size(test_trial_math_listen,2));
            for trial = 1:min_trial_count
                if sum(isnan(MSE_math_listen(trial,:)),2) == 0;
                [a max_count_math_listen(trial)] = nanmin(MSE_math_listen(trial,:)); % min SEM instead of MAX
                end
            end
            max_count_math_listen_correct = max_count_math_listen;
            max_count_math_listen_correct(max_count_math_listen_correct==2) = 0;
        
        % classification of task conditions with only math ANSWERING trials 
        % for each math test trial, calculate its correlation coefficient with  
        % 2 training templates (math listening and answering)
        for trial = 1:min_trial_count
            % remove voxels in both test trials and training templates if 0 
            % and nan values are found in either 
            % this step is to ensure the same data size with no nan/0 in 
            % all variables when calculating correlation coefficients
            temp1_test_math_answer = temp1_math_answer_test_angle(:,:,trial);
            temp1_train_math_listen_filt = temp1_train_math_listen(~isnan(temp1_test_math_answer));
            temp1_train_math_answer_filt = temp1_train_math_answer(~isnan(temp1_test_math_answer));    
            temp1_test_math_answer = temp1_test_math_answer(~isnan(temp1_test_math_answer));
            temp1_train_math_listen_filt_01 = temp1_train_math_listen_filt./temp1_train_math_listen_filt;
            temp1_train_math_answer_filt_01 = temp1_train_math_answer_filt./temp1_train_math_answer_filt;
            temp1_test_story_01 = temp1_test_math_answer./temp1_test_math_answer;    
            filt_notnan_01 =  temp1_train_math_listen_filt_01(:).*temp1_train_math_answer_filt_01(:).*temp1_test_story_01(:); 
            temp1_train_math_answer_filt_filt = temp1_train_math_answer_filt(~isnan(filt_notnan_01));
            temp1_train_math_listen_filt_filt = temp1_train_math_listen_filt(~isnan(filt_notnan_01));
            temp1_test_math_answer_filt = temp1_test_math_answer(~isnan(filt_notnan_01));

            if size(temp1_test_math_answer_filt(:)) == 1 % remove test trials with only 1 data points
                continue
            end
            % calcualte angle differences between phase vectors of
            % math answering test trials and math listening
            % training template
            error_answer_listen = temp1_test_math_answer_filt(:)-temp1_train_math_listen_filt_filt(:);
            % calcualte angle differences between phase vectors of
            % math answering test trials and math answering
            % training template     
            error_answer_answer = temp1_test_math_answer_filt(:)-temp1_train_math_answer_filt_filt(:);
            % restrict angle difference between -pi and pi
            for it = 1:size(error_answer_listen(:),1)
                if error_answer_listen(it)>pi
                    error_answer_listen(it) = error_answer_listen(it) - 2.*pi;
                elseif error_answer_listen(it) < -pi
                    error_answer_listen(it) = error_answer_listen(it) + 2.*pi;
                end     
            end
            for it = 1:size(error_answer_answer(:),1)
                if error_answer_answer(it)>pi
                    error_answer_answer(it) = error_answer_answer(it) - 2.*pi;
                elseif error_answer_answer(it) < -pi
                    error_answer_answer(it) = error_answer_answer(it) + 2.*pi;
                end     
            end            
            % calcualte mean-squared erros (MSE) of angle differences 
            % between phase vectors of math listening test trials
            % and math listening training template  
            MSE_answer_answer(trial) = nansum(error_answer_answer.^2)./nansum(error_answer_answer(:)./error_answer_answer(:));
            % calcualte mean-squared erros (MSE) of angle differences 
            % between phase vectors of math listening test trials
            % and math answering training template              
            MSE_answer_listen(trial) = nansum(error_answer_listen.^2)./nansum(error_answer_listen(:)./error_answer_listen(:));
        end
            % compare MSE between answer-answer and answer-listen pairs
            % (i.e., answer-listen = answering test trial vs. listeing training template)
            % the one with lower MSE is selected as the predicted task condition          
            MSE_math_answer = [MSE_answer_answer' MSE_answer_listen'];
            max_count_math_answer = nan(1,size(test_trial_math_answer,2));
            for trial = 1:min_trial_count
                if sum(isnan(MSE_math_answer(trial,:)),2) == 0;
                [a max_count_math_answer(trial)] = nanmin(MSE_math_answer(trial,:)); % min SEM instead of MAX
                end
            end
            max_count_math_answer_correct = max_count_math_answer;
            max_count_math_answer_correct(max_count_math_answer_correct==2) = 0;
            
            % combined prediction result of both math listen and answering test trials 
            % prediction accuracy as % of all predictions made, for each iteration 
            max_count_correct = [max_count_math_answer_correct max_count_math_listen_correct];
            max_count_correct_notnan = find(~isnan(max_count_correct)); % take out nan results
            max_count_correct_notnan_01 = max_count_correct_notnan./max_count_correct_notnan;
            classification_accuracy(iteration) = nansum(max_count_correct(:))./nansum(max_count_correct_notnan_01(:)); % 40 trials before, now depends on how many non-nan results

        end
        % iteraction-averaged prediction accuracy as % of all predictions 
        % made, for each voxel within the ROI
        classification_accuracy_avg(irow,icol) = nanmean(classification_accuracy(:));
        classification_accuracy_std(irow,icol) = nanstd(classification_accuracy(:));
        classification_accuracy_stderr(irow,icol) = classification_accuracy_std(irow,icol)./sqrt(iteration);
    end
end
%% save data
save_folder = [main_folder,'/Sample Data/Language Task Original 100 sub/Analysis/'];
save([save_folder,'classification_accuracy_local_phase_vector_field_MSE.mat'],'classification_accuracy_avg','classification_accuracy_stderr','classification_accuracy_std')   ;     
%% Visualization: 2D heatmap of classification accuracy of the local linear classifer
disp(['visualizing 2D heatmap of  classification accuracy...'])   
cd(main_folder)
load('parcellation_template7.mat')

figure()
pcolor(classification_accuracy_avg)
shading interp
colormap jet
colorbar
caxis([0.5,1])
ylim([1,176])
xlim([130,220])
title(['phase vector field angle based local (3x3) classification accuracy, math listening vs. answering'])
% overlayed with the borders of 7 parcellation template
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
hold off
disp(['completing process...'])    

end
