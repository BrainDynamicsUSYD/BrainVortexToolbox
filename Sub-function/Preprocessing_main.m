function [DataOut] = Preprocessing_main(subject,hemisphere,flagSur,flagRest,flagTask,flagSmooth)
disp(['initiating process...'])    
% add all the subfunctions
addpath(genpath([pwd]))
main_folder = pwd;
% set basic parameters
if flagSmooth == 0 && flagSur == 1 % no need to do surrogate raw data
    return
end
if flagSmooth == 2 && flagSur == 1 % no need to do smoothed surrogate data as it is done already when flagSmooth == 1
    return
end
%%
params.sigmScale =  [29.35 14.93];            % bandpass 5 bandwidth ranges
params.downSRate = 2 ;                        % downsample the re-interpolation
if hemisphere == 1                    
params.xCord = -250:params.downSRate:250 ;    % coordinate re-interpolation, left hemisphere
params.yCord = -150:params.downSRate:200 ;    
elseif hemisphere == 2                
params.xCord = -270:params.downSRate:230 ;    % coordinate re-interpolation, right hemisphere
params.yCord = -180:params.downSRate:170 ;    
end
params.fsTem = 1/0.72 ;                       % temporal sampling rate

% load raw and position data files
for iSub = subject:subject
    tic
    if flagRest == 1 % resting data
        name = dir([main_folder,'/Sample Data/Resting/Raw Data']) ;            % rest fMRI data  
        dataDir = [name(iSub+2).folder,'/',name(iSub+2).name,'/rfMRI_REST1_LR_Atlas.dtseries.nii'];    
        name = dir([main_folder,'/Sample Data/Resting/Data Pos']) ;        % position data of left/right hemisphere
        if hemisphere == 1  
            posFile = [name(iSub+2).folder,'/',name(iSub+2).name,'/','L.flat.32k_fs_LR.surf.gii'] ; % Left brain
        elseif hemisphere == 2  
            posFile = [name(iSub+2).folder,'/',name(iSub+2).name,'/','R.flat.32k_fs_LR.surf.gii'] ; % Right brain
        end               
    elseif flagRest == 0 % task data
        if flagTask == 1 % language task, original 100 subjects
        name = dir([main_folder,'/Sample Data/Language Task Original 100 sub/Raw Data']) ; % task data
        dataDir = [name(iSub+2).folder,'/',name(iSub+2).name,'/tfMRI_LANGUAGE_LR_Atlas.dtseries.nii']; 
        name = dir([main_folder,'/Sample Data/Language Task Original 100 sub/Data Pos']) ;        % position data of left/right hemisphere        
        if hemisphere == 1  
            posFile = [name(iSub+2).folder,'/',name(iSub+2).name,'/','L.flat.32k_fs_LR.surf.gii'] ; % Left brain
        elseif hemisphere == 2  
            posFile = [name(iSub+2).folder,'/',name(iSub+2).name,'/','R.flat.32k_fs_LR.surf.gii'] ; % Right brain
        end  
        elseif flagTask == 2 % language task, additional 100 subjects 
            name = dir([main_folder,'/Sample Data/Language Task Additional 100 sub/Raw Data']) ; % task data
            dataDir = [name(iSub+2).folder,'/',name(iSub+2).name,'/tfMRI_LANGUAGE_LR_Atlas.dtseries.nii']; 
            name = dir([main_folder,'/Sample Data/Language Task Additional 100 sub/Data Pos']) ;        % position data of left/right hemisphere        
            if hemisphere == 1  
                posFile = [name(iSub+2).folder,'/',name(iSub+2).name,'/',name(iSub+2).name,'.L.flat.32k_fs_LR.surf.gii'] ; % Left brain
            elseif hemisphere == 2  
                posFile = [name(iSub+2).folder,'/',name(iSub+2).name,'/',name(iSub+2).name,'.R.flat.32k_fs_LR.surf.gii'] ; % Right brain
            end  
        elseif flagTask == 3 % working memory task
            name = dir([main_folder,'/Sample Data/Working Memory Task/Raw Data']) ; % task data
            dataDir = [name(iSub+2).folder,'/',name(iSub+2).name,'/tfMRI_WM_LR_Atlas.dtseries.nii']; 
            name = dir([main_folder,'/Sample Data/Working Memory Task/Data Pos']) ;        % position data of left/right hemisphere        
            if hemisphere == 1  
                posFile = [name(iSub+2).folder,'/',name(iSub+2).name,'/','L.flat.32k_fs_LR.surf.gii'] ; % Left brain
            elseif hemisphere == 2  
                posFile = [name(iSub+2).folder,'/',name(iSub+2).name,'/','R.flat.32k_fs_LR.surf.gii'] ; % Right brain
            end      
        end
    end
    % process the data files
    cd(main_folder)
    surMethodNum = 7 ; flagVisBpSig = 0 ;   
    DataOut = load_fMRI(dataDir,posFile,flagSur,surMethodNum,params,flagVisBpSig,flagRest,flagSmooth,hemisphere,flagTask) ;

    % real data   
    if flagSur == 0
       if flagRest == 1
          if flagSmooth == 1 % temporal bandpass filtered
              if hemisphere == 1 % left hemisphere
                 folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                 cd(folder_name)
                 save(['Preprocessed_temporalbandpass_data_resting_LEFT_sub',num2str(iSub),'.mat'],'DataOut')   
              elseif hemisphere == 2 % right hemisphere
                 folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                 cd(folder_name)
                 save(['Preprocessed_temporalbandpass_data_resting_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')                 
              end
          elseif flagSmooth == 0 % raw data
              if hemisphere == 1 % left hemisphere
                 folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                 cd(folder_name)
                 save(['Preprocessed_raw_data_resting_LEFT_sub',num2str(iSub),'.mat'],'DataOut')   
              elseif hemisphere == 2 % right hemisphere
                 folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                 cd(folder_name)
                 save(['Preprocessed_raw_data_resting_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')                 
              end     
          elseif flagSmooth == 2 % spatiotemporal bandpss filtered data
              if hemisphere == 1 % left hemisphere
                 folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                 cd(folder_name)
                 save(['Preprocessed_spatiotemporalbandpass_data_resting_LEFT_sub',num2str(iSub),'.mat'],'DataOut')   
              elseif hemisphere == 2 % right hemisphere
                 folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                 cd(folder_name)
                 save(['Preprocessed_spatiotemporalbandpass_data_resting_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')                 
              end
          end
       elseif flagRest == 0
          if flagTask == 1 % language task, original 100 subjets
              if flagSmooth == 1 % temporal bandpass filtered
                  if hemisphere == 1 % left hemisphere
                     folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_temporalbandpass_data_language_task_orig100_LEFT_sub',num2str(iSub),'.mat'],'DataOut')   
                  elseif hemisphere == 2 % right hemisphere
                     folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_temporalbandpass_data_language_task_orig100_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')                 
                  end   
              elseif flagSmooth == 0 % raw data
                  if hemisphere == 1 % left hemisphere
                     folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_raw_data_language_task_orig100_LEFT_sub',num2str(iSub),'.mat'],'DataOut')   
                  elseif hemisphere == 2 % right hemisphere
                     folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_raw_data_language_task_orig100_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')                 
                  end  
              elseif flagSmooth == 2 % spatiotemporal bandpss filtered data
                  if hemisphere == 1 % left hemisphere
                     folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_spatiotemporalbandpass_data_language_task_orig100_LEFT_sub',num2str(iSub),'.mat'],'DataOut')   
                  elseif hemisphere == 2 % right hemisphere
                     folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_spatiotemporalbandpass_data_language_task_orig100_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')                 
                  end                   
              end
          elseif flagTask == 2 % language task, additional 100 subjets
              if flagSmooth == 1
                  if hemisphere == 1 % left hemisphere
                     folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_temporalbandpass_data_language_task_add100_LEFT_sub',num2str(iSub),'.mat'],'DataOut')   
                  elseif hemisphere == 2 % right hemisphere
                     folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_temporalbandpass_data_language_task_add100_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')                 
                  end          
              elseif flagSmooth == 0
                  if hemisphere == 1 % left hemisphere
                     folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_raw_data_language_task_add100_LEFT_sub',num2str(iSub),'.mat'],'DataOut')   
                  elseif hemisphere == 2 % right hemisphere
                     folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_raw_data_language_task_add100_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')                 
                  end    
              elseif flagSmooth == 2
                  if hemisphere == 1 % left hemisphere
                     folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_spatiotemporalbandpass_data_language_task_add100_LEFT_sub',num2str(iSub),'.mat'],'DataOut')   
                  elseif hemisphere == 2 % right hemisphere
                     folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_spatiotemporalbandpass_data_language_task_add100_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')                 
                  end                     
              end
          elseif flagTask == 3
              if flagSmooth == 1
                  if hemisphere == 1 % left hemisphere
                     folder_name = [main_folder,'/Sample Data/Working Memory Task/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_temporalbandpass_data_WM_task_LEFT_sub',num2str(iSub),'.mat'],'DataOut')   
                  elseif hemisphere == 2 % right hemisphere
                     folder_name = [main_folder,'/Sample Data/Working Memory Task/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_temporalbandpass_data_WM_task_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')                 
                  end                     
              elseif flagSmooth == 0
                  if hemisphere == 1 % left hemisphere
                     folder_name = [main_folder,'/Sample Data/Working Memory Task/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_raw_data_WM_task_LEFT_sub',num2str(iSub),'.mat'],'DataOut')   
                  elseif hemisphere == 2 % right hemisphere
                     folder_name = [main_folder,'/Sample Data/Working Memory Task/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_raw_data_WM_task_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')                 
                  end       
              elseif flagSmooth == 2
                  if hemisphere == 1 % left hemisphere
                     folder_name = [main_folder,'/Sample Data/Working Memory Task/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_spatiotemporalbandpass_data_WM_task_LEFT_sub',num2str(iSub),'.mat'],'DataOut')   
                  elseif hemisphere == 2 % right hemisphere
                     folder_name = [main_folder,'/Sample Data/Working Memory Task/Preprocessed Data'];
                     cd(folder_name)
                     save(['Preprocessed_spatiotemporalbandpass_data_WM_task_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')                 
                  end                      
              end
          end
           
       end
       
    % surrogate data (null model)  
    elseif flagSur == 1
       if flagRest == 1
          if flagSmooth == 1  % temporal and bandpass filtered (both smooth and unsmoothed surrogate data in one go to ensure consistency within the same batch of randomization)
              if hemisphere == 1 % left hemisphere
                 folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                 cd(folder_name)
                 DataOut_unsmooth = DataOut(:,:,:,2);  % unsmoothed
                 save(['Preprocessed_temporalbandpass_data_sur_resting_LEFT_sub',num2str(iSub),'.mat'],'DataOut_unsmooth')   
                 DataOut_smooth = DataOut(:,:,:,1);  % smoothed
                 save(['Preprocessed_spatiotemporalbandpass_data_sur_resting_LEFT_sub',num2str(iSub),'.mat'],'DataOut_smooth')   
              elseif hemisphere == 2 % right hemisphere
                 folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                 cd(folder_name)
                 DataOut_unsmooth = DataOut(:,:,:,2);  % unsmoothed
                 save(['Preprocessed_temporalbandpass_data_sur_resting_RIGHT_sub',num2str(iSub),'.mat'],'DataOut_unsmooth')   
                 DataOut_smooth = DataOut(:,:,:,1);  % smoothed
                 save(['Preprocessed_spatiotemporalbandpass_data_sur_resting_RIGHT_sub',num2str(iSub),'.mat'],'DataOut_smooth')                 
              end    
          end
       elseif flagRest == 0
          if flagTask == 1 % language task, original 100 subjets
              if flagSmooth == 1 % temporal and bandpass filtered
                  if hemisphere == 1 % left hemisphere
                     folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
                     cd(folder_name)
                     DataOut_unsmooth = DataOut(:,:,:,2);  % unsmoothed
                     save(['Preprocessed_temporalbandpass_data_sur_language_task_orig100_LEFT_sub',num2str(iSub),'.mat'],'DataOut_unsmooth')   
                     DataOut_smooth = DataOut(:,:,:,1);  % smoothed
                     save(['Preprocessed_spatiotemporalbandpass_data_sur_language_task_orig100_LEFT_sub',num2str(iSub),'.mat'],'DataOut_smooth')   
                  elseif hemisphere == 2 % right hemisphere
                     folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
                     cd(folder_name)
                     DataOut_unsmooth = DataOut(:,:,:,2);  % unsmoothed
                     save(['Preprocessed_temporalbandpass_data_sur_language_task_orig100_RIGHT_sub',num2str(iSub),'.mat'],'DataOut_unsmooth')   
                     DataOut_smooth = DataOut(:,:,:,1);  % smoothed
                     save(['Preprocessed_spatiotemporalbandpass_data_sur_language_task_orig100_RIGHT_sub',num2str(iSub),'.mat'],'DataOut_smooth')                  
                  end                   
              end
          elseif flagTask == 2 % language task, additional 100 subjets
              if flagSmooth == 1
                  if hemisphere == 1 % left hemisphere
                     folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Preprocessed Data'];
                     cd(folder_name)
                     DataOut_unsmooth = DataOut(:,:,:,2);  % unsmoothed
                     save(['Preprocessed_temporalbandpass_data_sur_language_task_add100_LEFT_sub',num2str(iSub),'.mat'],'DataOut_unsmooth')   
                     DataOut_smooth = DataOut(:,:,:,1);  % smoothed
                     save(['Preprocessed_spatiotemporalbandpass_data_sur_language_task_add100_LEFT_sub',num2str(iSub),'.mat'],'DataOut_smooth')                 
                  elseif hemisphere == 2 % right hemisphere
                     folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Preprocessed Data'];
                     cd(folder_name)
                     DataOut_unsmooth = DataOut(:,:,:,2);  % unsmoothed
                     save(['Preprocessed_temporalbandpass_data_sur_language_task_add100_RIGHT_sub',num2str(iSub),'.mat'],'DataOut_unsmooth')                 
                     DataOut_smooth = DataOut(:,:,:,1);  % smoothed
                     save(['Preprocessed_spatiotemporalbandpass_data_sur_language_task_add100_RIGHT_sub',num2str(iSub),'.mat'],'DataOut_smooth')                 
                  end                           
              end
          elseif flagTask == 3 % working memory task
              if flagSmooth == 1 
                  if hemisphere == 1 % left hemisphere
                     folder_name = [main_folder,'/Sample Data/Working Memory Task/Preprocessed Data'];
                     cd(folder_name)
                     DataOut_unsmooth = DataOut(:,:,:,2);  % unsmoothed
                     save(['Preprocessed_temporalbandpass_data_sur_WM_task_LEFT_sub',num2str(iSub),'.mat'],'DataOut_unsmooth')   
                     DataOut_smooth = DataOut(:,:,:,1);  % smoothed
                     save(['Preprocessed_spatiotemporalbandpass_data_sur_WM_task_RIGHT_sub',num2str(iSub),'.mat'],'DataOut_smooth')                 
                  elseif hemisphere == 2 % right hemisphere
                     folder_name = [main_folder,'/Sample Data/Working Memory Task/Preprocessed Data'];
                     cd(folder_name)
                     DataOut_unsmooth = DataOut(:,:,:,2);  % unsmoothed
                     save(['Preprocessed_temporalbandpass_data_sur_WM_task_RIGHT_sub',num2str(iSub),'.mat'],'DataOut_unsmooth')                 
                     DataOut_smooth = DataOut(:,:,:,1);  % smoothed
                     save(['Preprocessed_spatiotemporalbandpass_data_sur_WM_task_RIGHT_sub',num2str(iSub),'.mat'],'DataOut_smooth')                 
                  end                                              
              end
          end
           
       end       
        
        
    end
    


        disp(['finishing subject ',num2str(iSub)])
        toc
        cd(main_folder)
    end
end

