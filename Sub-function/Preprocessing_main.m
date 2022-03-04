function [sigBPass] = Preprocessing_main(subject,hemisphere,flagSur)
disp(['initiating process...'])    
% add all the subfunctions
addpath(genpath([pwd]))

% set basic parameters
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
params.zscore = 0 ;                           % flag, 1 for zscore, 0 for no     

    for iSub = subject:subject
        tic
        name = dir([pwd,'/Sample Data/Raw Data']) ;            % task fMRI data  
        dataDir = [name(iSub+2).folder,'/',name(iSub+2).name,'/tfMRI_LANGUAGE_LR_Atlas.dtseries.nii']; 

        name = dir([pwd,'/Sample Data/Data Pos']) ;        % position data of left/right hemisphere
        if hemisphere == 1  
        posFile = [name(iSub+2).folder,'/',name(iSub+2).name,'/','L.flat.32k_fs_LR.surf.gii'] ; % Left brain
        elseif hemisphere == 2  
        posFile = [name(iSub+2).folder,'/',name(iSub+2).name,'/','R.flat.32k_fs_LR.surf.gii'] ; % Right brain
        end
        
        surMethodNum = 3 ; flagVisBpSig = 0 ; % flagSur = 1 to produce surrogate data with surrogate method 'surMethodNum = 1-6' 
        sigBPass = load_fMRI(dataDir,posFile,flagSur,surMethodNum,params,flagVisBpSig) ;
        save_folder = [pwd,'/Sample Data/'];
        save([save_folder,'Preprocessed_bandpass_ScaleSigma_29.35_14.93_sub',num2str(iSub),'.mat'],'sigBPass')        
        disp(['finishing subject ',num2str(iSub)])
        toc
    end
end

