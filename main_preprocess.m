% add all the subfunctions (folder preprocess)
% addpath(genpath([pwd,'/cifti-matlab-master'])) % cifti to read data !!!!
addpath(genpath([pwd,'/sub_fun']))

% set basic parameters
% params.sigmScale = 1.7 ; params.flagBP = 0;             % FHWM/2.35/2= sigma; 8/2.35/2 = 1.7
% params.sigmScale = 2.^(6:-1:2) ; params.flagBP = 2;
% params.sigmScale = [83.68,42.99,29.35,22.23,17.87,14.93] ;%,12.82,11.22,9.98,8.99,8.17,7.49] ; params.flagBP = 1;
% params.sigmScale = [83.68,29.35,17.87,12.82,9.98,8.17] ; params.flagBP = 1;
params.sigmScale = [29.35,14.93] ;params.flagBP = 1;
params.downSRate = 2 ;                        % downsample the re-interpolation
params.xCord = -250:params.downSRate:250 ;    % coordinate re-interpolation
params.yCord = -150:params.downSRate:200 ;
% params.xCord = -270:params.downSRate:230 ;    % coordinate re-interpolation
% params.yCord = -180:params.downSRate:170 ;
params.fsTem = 1/0.72 ;                       % temporal sampling rate
params.zscore = 0 ;                           % flag, 1 for zscore, 0 for no     
params.flagSur =1  ; params.surMethodNum =6 ; 
params.flagVisBpSig =0 ; 

% load([pwd,'/Results_data/fMRI_data/maskInfo.mat'])
for iSub = 1%:100%:length(name)-2
    % name = dir([pwd,'/Data/fMRI_Converted']) ;
    % name = dir([pwd,'/Data/fMRI_converted_motor']) ;
    name = dir([pwd,'/Data/fMRI_converted_language']) ;
    tic
    % dataDir = [name(iSub+2).folder,'/',name(iSub+2).name,'/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii'] ;
    % dataDir = [name(iSub+2).folder,'/',name(iSub+2).name,'/tfMRI_MOTOR_LR_Atlas_MSMAll.dtseries.nii'] ;
    dataDir = [name(iSub+2).folder,'/',name(iSub+2).name,'/tfMRI_LANGUAGE_LR_Atlas.dtseries.nii'] ;
    name = dir([pwd,'/Data/fMRI_Converted']) ;
    posFile = [name(iSub+2).folder,'/',name(iSub+2).name,'/L.flat.32k_fs_LR.surf.gii'] ;
    % posFile = ['/import/headnode2/xlon3884/Data/fMRI2/100610_3T_structural/MNINonLinear/fsaverage_LR32k/100610.R.flat.32k_fs_LR.surf.gii'] ;
    % flagSur = 1 ; surMethodNum = 4 ; flagVisBpSig = 0 ;
    params.name = name(iSub+2).name ;
    sigBPass = load_fMRI_pp(dataDir,posFile,params) ; 
    % save([pwd,'/preprocess/Preprocess_save/Preprocessed_sub',num2str(iSub),'.mat'],'sigBPass','params')
    % save([pwd,'/preprocess/Preprocess_LoG/Preprocessed_sub',num2str(iSub),'.mat'],'sigBPass','params')
    %save([pwd,'/preprocess/Preprocess_eigen345_right/Preprocessed_sub',num2str(iSub),'.mat'],'sigBPass','params')
    % save([pwd,'/preprocess/Preprocess_lp_zscore/Preprocessed_sub',num2str(iSub),'.mat'],'sigBPass','params')
    % save([pwd,'/preprocess/Preprocess_lp_1.7/Preprocessed_sub',num2str(iSub),'.mat'],'sigBPass','params')
    % save([pwd,'/preprocess/Preprocess_surBP_3spaceTime/Preprocessed_sub',num2str(iSub),'.mat'],'sigBPass','params')
    % save([pwd,'/preprocess/Preprocess_surrogate/Prepro_sur_sub',num2str(iSub),'.mat'],'sigBPass','params')
    % save([pwd,'/preprocess/Preprocess_motor/Preprocess_motor',num2str(iSub),'.mat'],'sigBPass','params')
    % save([pwd,'/preprocess/Preprocess_surBP/Prepro_sur_sub',num2str(iSub),'.mat'],'sigBPass','params')
    % save([pwd,'/preprocess/Preprocess_eigen345_lan/Preprocessed_sub',num2str(iSub),'.mat'],'sigBPass','params')
    disp(['finishing subject ',num2str(iSub)])
    toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% extract the vortex information
% testVorRegBen
cd preprocess
paramsVor.curlThre = 1 ;
paramsVor.thresholdTime = 10 ;
paramsVor.thresholdSize = 3 ;
paramsVor.expanPara = 1 ;

% fileDir = [pwd,'/Preprocess_save/'] 
fileDir = [pwd,'/Preprocess_eigen345_lan/'] ;
% fileDir = [pwd,'/Preprocess_lp_1.7/'] ;
% fileDir = [pwd,'/Preprocess_surrogate/'] ;
% fileDir = [pwd,'/Preprocess_surBP/'] ;
% fileDir = [pwd,'/Preprocess_BP/'] ;
bandNum = 1 ;
saveFolderName = 'VortexStat_eigen345_lan' ;


fileName1 = dir(fileDir) ;
for iSub = 1:length(fileName1)-2 ;
    tic
fileNum = 2 + iSub ;  % find the resting data

tic
load([fileDir,fileName1(fileNum).name])
sigIn = squeeze(sigBPass(bandNum,:,:,:)) ;
% sigInReshape = reshape(sigIn,176*251,[]) ;
% temp = find(~isnan(sigInReshape(:,1))) ;
% randSeq = randperm(length(temp)) ;
% sigInReshape(temp(randSeq),:) = sigInReshape(temp,:) ;
% sigIn = reshape(sigInReshape,176,251,[]) ;

hilbertSig = [] ;
phaseSig = [] ;

for iX = 1:size(sigIn,1)
    for iY = 1:size(sigIn,2)
        hilbertSig(iX,iY,:) = hilbert(sigIn(iX,iY,:)) ;
        phaseSig(iX,iY,:) = angle(hilbertSig(iX,iY,:)) ;
    end
end
toc

% find the phase gradient
tic
vPhaX = zeros(size(phaseSig)) ;
vPhaY = zeros(size(phaseSig)) ;
for iTime2 = 1:size(phaseSig,3)
    for iX = 1:size(phaseSig,1)
        vPhaX(iX,2:end-1,iTime2) = anglesubtract(phaseSig(iX,3:end,iTime2),phaseSig(iX,1:end-2,iTime2))/2 ;
    end
    for iY = 1:size(phaseSig,2)
        vPhaY(2:end-1,iY,iTime2) = anglesubtract(phaseSig(3:end,iY,iTime2),phaseSig(1:end-2,iY,iTime2))/2 ;
    end
end
toc
temp = vPhaX ;
vPhaseX = -vPhaX./sqrt(vPhaX.^2+vPhaY.^2) ;
vPhaseY = -vPhaY./sqrt(temp.^2+vPhaY.^2) ;
% Vx_flowmap_norm_phase = vPhaseX ;
% Vy_flowmap_norm_phase = vPhaseY ;
toc

%% find the vortex regions
tic
disp('start detecting vortex')
curlz = [] ;
cav = [] ;
for time = 1:size(sigIn,3)
    temp1_vx = vPhaseX(:,:,time);
    temp1_vy = vPhaseY(:,:,time);
    [x,y] = meshgrid(1:size(temp1_vx,2),1:size(temp1_vx,1));
    
    [curlz(:,:,time) cav(:,:,time)]  = curl(x,y,temp1_vx,temp1_vy);
end
cav_filt_pos = curlz;
cav_filt_pos(cav_filt_pos<paramsVor.curlThre) = 0;
%         cav_filt(cav_filt>0.6) = 1;
cav_filt_pos(isnan(cav_filt_pos)) = 0;

cav_filt_neg = curlz;
cav_filt_neg(cav_filt_neg>-paramsVor.curlThre) = 0;
%         cav_filt(cav_filt>0.6) = 1;
cav_filt_neg(isnan(cav_filt_neg)) = 0;

[durationValid,~,WCenPos,pattStPos]= detectCoM(cav_filt_pos,curlz,paramsVor) ;

% 2 for positive vortex and 1 for negative vortex
vortex_filt_pos = findVortexRegion(durationValid,WCenPos,pattStPos,vPhaseX,vPhaseY,curlz,paramsVor) ;

[durationValid,~,WCenNeg,pattStNeg]= detectCoM(cav_filt_neg,curlz,paramsVor) ;

vortex_filt_neg = findVortexRegion(durationValid,WCenNeg,pattStNeg,vPhaseX,vPhaseY,curlz,paramsVor) ;

%% combine negative and positive
vortex_template_posi_nega_sparse = [] ;
tic
% for subject = 1:No_subjects
%     subject
    for time = 1:max(size(vortex_filt_pos,2),size(vortex_filt_neg,2))
        
        vortex_template = nan(size(cav_filt_neg,1),size(cav_filt_neg,2),2);
        
        % positve only vortex
        if time<=size(vortex_filt_pos,2)
        for ipatt = 1:size(vortex_filt_pos,1);
            temp1 = vortex_filt_pos{ipatt,time};
            if nansum(temp1(:)) == 0;
                continue
            else
                temp_full = full(temp1);
                [row col] = find(temp_full);
                
                for i = 1:size(row,1)
                    vortex_template(row(i),col(i),1) = temp_full(row(i),col(i)) ;
                end
            end
        end
        end
        % negative only vortex
        if time<=size(vortex_filt_neg,2)
        for ipatt = 1:size(vortex_filt_neg,1);
            temp1 = vortex_filt_neg{ipatt,time};
            if nansum(temp1(:)) == 0;
                continue
            else
                temp_full = full(temp1);
                [row col] = find(temp_full);
                
                for i = 1:size(row,1)
                    vortex_template(row(i),col(i),2) = temp_full(row(i),col(i)) ;
                end
            end
        end
        end
        % vortex_template(vortex_template==2) = 1; % to cancel out pos and
        % neg
        vortex_template_posi_nega = reshape(nansum(vortex_template,3),[size(cav_filt_neg,1),size(cav_filt_neg,2)]);
        vortex_template_posi_nega(isnan(vortex_template_posi_nega)) = 0;
        vortex_template_posi_nega_sparse{time} = vortex_template_posi_nega;
        
    end
% end

toc

% % visualization
%     sigSmooth = [] ;
%     resizeScale = 1 ;
%     sigIn3 = curlz ;
%     for iTime2 = 1:size(sigIn3,3)
%         %Try smoothing
%         filtWidth = 5;
%         filtSigma = 0.6;   % 0.6
%         imageFilter=fspecial('gaussian',filtWidth,filtSigma);
%         smoothTemp = nanconv((squeeze(sigIn3(:,:,iTime2))),imageFilter,'edge', 'nonanout');
%         sigSmooth(:,:,iTime2) = imresize(smoothTemp, resizeScale);
%     end

%% visualize vortex full
flagPlot = 0 ;
if flagPlot
flagSaveVideo = 0 ;
if flagSaveVideo
    t =  datetime('now') ;
    dateStr = datestr(t,'mmmmdd_HH:MM') ;
    vidTitle = [pwd,'/Results/Project3/movies/'] ;
    saveFileName = [vidTitle,'fMRI_s1','vortex_regionVortSm',dateStr] ;
    vidObj = VideoWriter(saveFileName,'Motion JPEG AVI');
    v.Quality = 50 ;
    vidObj.FrameRate = 2 ;  % 5
    open(vidObj);
    fig = figure;
    set(gcf,'Position', [360 156 710 530])
    % set(gcf,'Position',[358 266 560 420])
    % set(gcf,'Position',[39 251 1168 435]) % for 2 subplot
else
    figure;
end
m = 256 ;

for iTime = 1:1200
    % subplot(1,2,1)
    binPatt = vortex_template_posi_nega_sparse{iTime} ;
    binPatt(binPatt~=0) = 1 ;
    imagesc(binPatt.*(sigSmooth(:,:,iTime)))
    % imagesc(vortex_template_posi_nega_sparse{iTime}(:,:))
    caxis([-1,1])
    % c = [0,0,1;1,1,1;1,1,0;1,0,0] ;
        m1 = m*0.25;
    r = [zeros(m1,1); (0:m1-1)'/max(m1-1,1)] ;
    r2 = [ones(m1,1); (m1-1:-1:0)'/max(m1-1,1)/2+0.5 ] ;
    r = [r; r2];
    g = [zeros(m1,1);(0:m1-1)'/max(m1-1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
    
c = [r g b];
% c(1,:) = 1 ;
    colormap(c)
    set(gca,'YDir','normal')
    downSR = 3 ;
    hold on
    quiver(1:downSR:251,1:downSR:176,vPhaseX(1:downSR:end,1:downSR:end,iTime)...
        ,vPhaseY(1:downSR:end,1:downSR:end,iTime)...
        ,0.4 ,'Color',[0.4 0.4 0.4], 'MaxHeadSize',1)
    
    plot22SecBound
    load([pwd,'/Results_data/fMRI_data/maskInfo.mat'])
    b_plot = [(b(:,1)-min(params.xCord)+1)/2+0.5, (b(:,2)-min(params.yCord)+1)/2+0.5 ] ;
    lh = plot(b_plot(:,1),b_plot(:,2),'k--','linewidth',0.5) ;
    lh.Color = [lh.Color 0.8] ;
    title(['vortices at time frame = ',num2str(iTime)])
    
%     visX = 75:125 ;  visY = 140:200 ;
%     subplot(1,2,2)
%     % imagesc(visY,visX,vortex_template_posi_nega_sparse{iTime}(visX,visY))
%     imagesc(visY,visX,binPatt(visX,visY).*(sigSmooth(visX,visY,iTime)))
%     set(gca,'YDir','normal')
%     caxis([-1,1])
%     % c = [0,0,1;1,1,1;1,1,0;1,0,0] ;
%     colormap(c)
%     hold on
%     quiver(visY,visX,vPhaseX(visX,visY,iTime),vPhaseY(visX,visY,iTime)...
%         ,0.4 ,'Color',[0.4 0.4 0.4], 'MaxHeadSize',1)
%     
%     plot22SecBound
% %     load([pwd,'/Results_data/fMRI_data/maskInfo.mat'])
% %     b_plot = [(b(:,1)-min(params.xCord)+1)/2+0.5, (b(:,2)-min(params.yCord)+1)/2+0.5 ] ;
% %     lh = plot(b_plot(:,1),b_plot(:,2),'k--','linewidth',0.5)
% %     lh.Color = [lh.Color 0.8] ; 
%     
%     title(['LP at time frame = ',num2str(iTime)])
    if flagSaveVideo
        writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
    else
        pause
    end
    cla
end
if flagSaveVideo
    close(vidObj)
end

end

%% study the spatial distribtuion
newVor = zeros(size(cav_filt_neg,1),size(cav_filt_neg,2),size(vortex_filt_pos,2)) ;
for iTime = 1:size(vortex_template_posi_nega_sparse,2)
    newVor(:,:,iTime) = vortex_template_posi_nega_sparse{iTime} ;
end
flagSave = 1 ;
if flagSave
% save([pwd,'/Results_data/fMRI_data/vorRegion_lan/sub',num2str(iSub),'.mat'],...
%     'newVor','pattStPos','pattStNeg','WCenNeg','WCenPos','vortex_filt_neg','vortex_filt_pos','paramsVor') ;
% save([pwd,'/VortexStat_save/Vortex_sub',num2str(iSub),'.mat'],...
%     'newVor','pattStPos','pattStNeg','WCenNeg','WCenPos','vortex_filt_neg','vortex_filt_pos','paramsVor') ;
% save([pwd,'/VortexStat_LoG3/Vortex_sub',num2str(iSub),'.mat'],...
%     'newVor','pattStPos','pattStNeg','WCenNeg','WCenPos','vortex_filt_neg','vortex_filt_pos','paramsVor') ;
% save([pwd,'/VortexStat_lp_1.7/Vortex_sub',num2str(iSub),'.mat'],...
%     'newVor','pattStPos','pattStNeg','WCenNeg','WCenPos','vortex_filt_neg','vortex_filt_pos','paramsVor') ;
% save([pwd,'/vortexStat_sur/Vortex_sub',num2str(iSub),'.mat'],...
%     'newVor','pattStPos','pattStNeg','WCenNeg','WCenPos','vortex_filt_neg','vortex_filt_pos','paramsVor') ;
% save([pwd,'/vortexStat_surBP_band4/Vortex_sub',num2str(iSub),'.mat'],...
%    'newVor','pattStPos','pattStNeg','WCenNeg','WCenPos','vortex_filt_neg','vortex_filt_pos','paramsVor') ;
% save([pwd,'/VortexStat_band2/Vortex_sub',num2str(iSub),'.mat'],...
%     'newVor','pattStPos','pattStNeg','WCenNeg','WCenPos','vortex_filt_neg','vortex_filt_pos','paramsVor') ;
save([pwd,'/',saveFolderName,'/Vortex_sub',num2str(iSub),'.mat'],...
    'newVor','pattStPos','pattStNeg','WCenNeg','WCenPos','vortex_filt_neg','vortex_filt_pos','paramsVor') ;

end
newVor(newVor~=0) = 1 ;
spaceDist = nansum(newVor,3) ;
spaceDistAll(:,:,iSub) = spaceDist ;
disp(['finishing subject ',num2str(iSub)])
toc

end

