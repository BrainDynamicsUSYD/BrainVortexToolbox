%% for Manuscript figures;
% Xian Long  xlon3884@uni.sydney.edu.au
% 08/08/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% global variables
% testVorRegBen
addpath(genpath([pwd,'/Toolbox_CSC']))
addpath(genpath([pwd,'/ToolOthers/nanconv']))
addpath(genpath([pwd,'/ToolOthers/fMRI/cifti-matlab-master']))
addpath(genpath([pwd,'/ToolOthers/fMRI/BrainSpace-master/matlab']))

% set basic parameters
params.sigmScale = 2.^(6:-1:2) ;              % log spacing of scales
params.downSRate = 2 ;                        % downsample the re-interpolation
params.xCord = -250:params.downSRate:250 ;    % coordinate re-interpolation
params.yCord = -150:params.downSRate:200 ;
params.fsTem = 1/0.72 ;                       % temporal sampling rate
params.zscore = 1 ;  

%% preprocessing 
iSub = 1 ;
% fileDir = ['/Results_data/fMRI_data/Preprocess_large/'] ;
% fileName1 = dir([pwd,'/Data/fMRI_Converted/']) ;
% fileNum = 2 + iSub ;  % find the resting data
% fileName2 = dir([pwd,fileDir,'*rest*',fileName1(fileNum).name,'*']) ;
% % load pre-processed data
% load([fileName2.folder,'/',fileName2.name])
load([pwd,'/preprocess/Preprocess_eigen345/Preprocessed_sub',num2str(iSub),'.mat'])
sigIn = squeeze(sigBPass) ;

% % load vortex dynamics
% load([pwd,'/Results_data/fMRI_data/vorRegion_rest_0.5expan/sub',num2str(iSub)]) ;
load([pwd,'/preprocess/VortexStat_eigen345/Vortex_sub',num2str(iSub),'.mat'])

flagFindPhase = 1 ;
if flagFindPhase
    tic
    hilbertSig = [] ;
    phaseSig = [] ;
    
    for iX = 1:size(sigIn,1)
        for iY = 1:size(sigIn,2)
            hilbertSig(iX,iY,:) = hilbert(sigIn(iX,iY,:)) ;
            phaseSig(iX,iY,:) = angle(hilbertSig(iX,iY,:)) ;
        end
    end
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
    temp = vPhaX ;
    vPhaseX = -vPhaX./sqrt(vPhaX.^2+vPhaY.^2) ;
    vPhaseY = -vPhaY./sqrt(temp.^2+vPhaY.^2) ;
    toc
end

%% fig. 2A. 3D brain and time series (amp. & phase)
% further bandpass for visualization
fLow = 0.01 ;
fHigh = 0.05 ;
sigInReshape = reshape(sigIn,176*251,[]) ;
[bandpasSigRS,~, ~, ~] = bandpa_fMRI(sigInReshape,params.fsTem ,fLow,fHigh) ;    
bandpasSig = reshape(bandpasSigRS,176,251,1200) ;

phaseSigBP_RS = angle(hilbert(bandpasSigRS'))' ;
phaseSigBP = reshape(phaseSigBP_RS,176,251,1200) ;


%%
% close all
figure;
% timeSec = 600:800 ;
timeSec = 240:440 ;
% timeSec = 620:660 ;
x1 = 88 ; y1 = 182 ;x2 = 98 ; y2 = 162 ;  % Opposite
% x1 = 88 ; y1 = 182 ;x2 = 88 ; y2 = 172 ;  % lagging

plot(timeSec,squeeze(bandpasSig(x1,y1,timeSec)),'b','lineWidth',2)
xlim([timeSec(1) timeSec(end)])
%ylim([-0.8 0.8])
hold on
plot(timeSec,squeeze(bandpasSig(x2,y2,timeSec)),'r','lineWidth',2)
xlim([timeSec(1) timeSec(end)])
% ylim([-0.6 0.6])

figure;
plot(timeSec,squeeze(phaseSigBP(x1,y1,timeSec)),'b','lineWidth',2)
xlim([timeSec(1) timeSec(end)])
ylim([-pi pi])
hold on
plot(timeSec,squeeze(phaseSigBP(x2,y2,timeSec)),'r','linewidth',2)
xlim([timeSec(1) timeSec(end)])
ylim([-pi pi])

%% whole brain
folderName = [pwd,'/Data/fMRI2/100610_3T_structural/MNINonLinear/fsaverage_LR32k/'] ;
fileName = '100610.L.inflated.32k_fs_LR.surf.gii' ;
posFile = [folderName,fileName] ;
posLC_inflat = gifti(posFile) ;

    viewAngle = [-90 0] ;
    temp = ones(32492,1) ;
    clim = [0,1] ;
    h = [] ;
    figure;
    convert2Inflated(temp,posLC_inflat,[],clim,viewAngle)
    colormap(ones(256,1)*[1,0.6,0.6])
% hold on
% plot3(posLC_inflat.vertices(idx,1),posLC_inflat.vertices(idx,2)...
%          ,posLC_inflat.vertices(idx,3),'--','Color',[0.4,0.4,0.4])
h.camplight= camlight() ;


%% fig. 2BC. patt whole brain
% close all

visX = 1:176 ; visY = 1:251 ;
vDownR = 3 ;
saveName = 'wholeBrain' ;
genFig1A(sigIn,visX,visY,vDownR,0,saveName) ;

% fig. 1B. patt lateral parietal
%visX = 84:108 ; visY = 136:160 ;
visX = 96:124 ; visY = 168:196 ;

vDownR = 1 ;
saveName = 'regionLP' ;
genFig1A(sigIn,visX,visY,vDownR,0,saveName) ;

%% fig. 3A. propagating vortices
close all
for iSub = 1%:60
    newVor = zeros(176,251,1200) ;
    % load([pwd,'/Results_data/fMRI_data/VorRegionFull/vorRegion_',num2str(iSub),'.mat']) ;
    % load([pwd,'/Results_data/fMRI_data/VorRegionFull/vorRegion_rest/sub',num2str(iSub),'.mat']) ;
    pattStPos = pattStPos ;
wCenPos = WCenPos ;
flagStepOne = 1 ;
    for iCentPos = 73 %90% :size(WCenPos,2)
        iPlot = 1 ;
        plotSt =  4% 6 ;
        plotGap =  4 ; %3
        plotEnd = plotSt+plotGap*3 ;
        if size(WCenPos{iCentPos},1)<plotEnd
            continue
        end
        figure
        for iTime = pattStPos(iCentPos)+plotSt-1:plotGap:pattStPos(iCentPos)+plotEnd-1 %12%21%28 pattStPos(iCentPos):pattStPos(iCentPos)+size(WCenPos{iCentPos},1)-1
            subplot(1,4,iPlot);
            iPlot = iPlot + 1 ;
            curReg = vortex_filt_pos{iCentPos,iTime} ;
            curReg(curReg==0) = nan ;
            [I,J] = find(~isnan(curReg)) ;
            radSet = ceil(min([max(J)-min(J),max(I)-min(I)])/2 ) ;
            [X,Y] = meshgrid(1:radSet*2,1:radSet*2);
            c=[radSet+0.5 radSet+0.5];%center coordinates, these are allowed to be on the edge or even outside of the image
            dis_loc = double(sqrt((X-c(1)).^2+(Y-c(2)).^2) <= radSet);
            dis_loc(dis_loc==0) = nan ;
            
            curCen = iTime-pattStPos(iCentPos) ;  % weighted center time

            curReg = nan(size(curReg)) ;
            curReg(ceil(WCenPos{iCentPos}(curCen,2))-(radSet):floor(WCenPos{iCentPos}(curCen,2))+(radSet),...
                ceil(WCenPos{iCentPos}(curCen,1))-(radSet):floor(WCenPos{iCentPos}(curCen,1))+(radSet))...
                = dis_loc ;
            curPhase = curReg.*phaseSig(:,:,iTime) ;
            imagesc(curPhase)
            colormap([1,1,1;parula(255)])
            % axis off
            % xlim([90 130]) ; ylim([140 170])  % patt = 1
            % xlim([110 150]) ; ylim([60 100])  % patt = 2
            % xlim([120 150]) ; ylim([130 160])  % patt = 3
            % xlim([130 160]) ; ylim([40 70])  % patt = 4
            % xlim([170 210]) ; ylim([30 60])  % patt = 5
            % xlim([184 224]) ; ylim([50 90])  % patt = 6
            % xlim([170 210]) ; ylim([120 160])  % patt = 7
            % xlim([120 150]) ; ylim([100 130])  % patt = 14
            % xlim([120 150]) ; ylim([130 160])  % patt = 22
            % xlim([210 250]) ; ylim([70 100])  % patt = 26

            hold on;
            plot(WCenPos{iCentPos}(1:curCen,1),WCenPos{iCentPos}(1:curCen,2),'r','lineWidth',1)
            plot(WCenPos{iCentPos}(curCen,1),WCenPos{iCentPos}(curCen,2),'r.','markerSize',12)
            % xlim([WCenPos{iCentPos}(1,1)-6 WCenPos{iCentPos}(1,1)+14]) ; 
            % ylim([WCenPos{iCentPos}(1,2)-4 WCenPos{iCentPos}(1,2)+16])
            xlim([WCenPos{iCentPos}(1,1)-10 WCenPos{iCentPos}(1,1)+12]) ; 
            ylim([WCenPos{iCentPos}(1,2)-10 WCenPos{iCentPos}(1,2)+12])
            caxis([-pi,pi])
            set(gca,'yDir','normal')
            % pause
            % cla
            curPhase = curReg.*phaseSig(:,:,iTime) ;
            [I,J] = find(abs(curPhase-pi)<0.2) ;      % find the coord. wh phase equal pi
            iTime2 = iTime - pattStPos(iCentPos) ;  % for centres
            % plot(J,I,'r.')    % check
            if isempty(I)
                continue
            end
            if flagStepOne
                flagStepOne = 0 ;
                angStepOne = mean((angle(J - WCenPos{iCentPos}(iTime2,1)+1+ ...
                    i*(I - WCenPos{iCentPos}(iTime2,2)+1 )))) ;
            end
            curRad = mean((angle(J - WCenPos{iCentPos}(iTime2,1)+ ...
                    i*(I - WCenPos{iCentPos}(iTime2,2))))) ;
            AngVel = (mod( (curRad)-angStepOne ,2*pi)) ; 
            if AngVel>5
                AngVel = 0 ;
            end
            title(['\theta = ',num2str(round(rad2deg(curRad)),'%d')...
                ', \omega = ',num2str(round(rad2deg(AngVel)),'%d')])
            
        end
    end

end

% exportgraphics(gcf,['3A.pdf'],'ContentType','vector')

%% fig. 3B. many vortices


%% fig. 3C. MSD example
totalValid = [] ;
totalSuper = [] ;
totalBurst = [] ;

for iSub = 1:180
    load([pwd,'/preprocess/VortexStat_eigen345/Vortex_sub',num2str(iSub),'.mat'])
    [~,meanSlope,stdSlope,slopeAll,numPatt] = fMRI_msd(WCenPos) ;
    meanSlope
    stdSlope
    totalValid(iSub) = length(slopeAll) ;
    totalSuper(iSub) = length(find(slopeAll>1)) ;
    totalBurst(iSub) = numPatt ;
    % fMRI_msd_example
    [~,meanSlope,~,slopeAll] = fMRI_msd(WCenNeg) ;
    meanSlope
    stdSlope
    totalValid(iSub+180) = length(slopeAll) ;
    totalSuper(iSub+180) = length(find(slopeAll>1)) ;
    totalBurst(iSub+180) = numPatt ;
        disp(['finishing subj. ',num2str(iSub)])

end


%% fig. 4A-E. vortex interaction
% required pre-processed raw data, vortex trajectories.
% function 'figInteract'
iSub = 1 ;
% disp('loading pre-processed data ...')
tic
fileDir = ['/Results_data/fMRI_data/Preprocess_large/'] ;
fileName1 = dir([pwd,'/Data/fMRI_Converted/']) ;
fileNum = 2 + iSub ;  % find the resting data
fileName2 = dir([pwd,fileDir,'*rest*',fileName1(fileNum).name,'*']) ;
load([fileName2.folder,'/',fileName2.name])
load([pwd,'/Results_data/fMRI_data/vorRegion_rest/sub',num2str(iSub),'.mat']) ;
toc

stdPara = 3 ;
% A annilihalation
%figInteract(sigIn,stdPara,WCenPos,pattStPos,WCenNeg,pattStNeg,params,1,'DipoleAnni') ;
% B partial annilihalation
figInteract(sigIn,stdPara,WCenPos,pattStPos,WCenNeg,pattStNeg,params,2,'DipolePartial') ;
% C repulsion dipole
figInteract(sigIn,stdPara,WCenPos,pattStPos,WCenNeg,pattStNeg,params,3,'DipoleRepul') ;
% D chain interaction
%figInteract(sigIn,stdPara,WCenPos,pattStPos,WCenNeg,pattStNeg,params,4,'chain') ;
% E repulsion same sign
figInteract(sigIn,stdPara,WCenPos,pattStPos,WCenNeg,pattStNeg,params,5,'SameRepul') ;

%% fig. 5A formation mechanism of vortices 1
% close all
%for iTime = 20:10:200
stepSize = 3 ;
%visX2 = 72:116 ;  visY2 = 148:192 ;   iTime = 262 ;   % LP
% visX2 = 64:116 ;  visY2 = 140:192 ; iTime = 628 ; 
visX2 =  128:156 ;  visY2 = 178:202 ;  iTime = 272 ;  % PCC  % visX2 =  110:158 ;  visY2 = 172:220 ;
flagSavePdf = 1 ;

for iTime = [270,273,274,275,278] ;
explainVor(sigIn,vPhaseX,vPhaseY,newVor,visX2,visY2,iTime,stepSize,flagSavePdf)
end
% explainVorSche(sigIn,vPhaseX,vPhaseY,newVor,visX2,visY2,iTime,stepSize,flagSavePdf)
% end

%% fig. 5C. spatial distribution of vortices
spaDist = zeros(176,251) ;          % centre dist
spaDistlong = zeros(176,251) ; 
spaVorDist = zeros(176,251) ;       % vorticities dist
spaFullVorDist = zeros(176,251) ;  % full vortex dist

for iSub = 1:180
tic
    load([pwd,'/preprocess/VortexStat_eigen345New/Vortex_sub',num2str(iSub),'.mat'])
    % load([pwd,'/preprocess/VortexStat_eigen345_rightNew/Vortex_sub',num2str(iSub),'.mat'])
    newVor(newVor~=0) = 1;
    spaFullVorDist = spaFullVorDist+sum(newVor,3) ;  % full vortex dist

    for iPatt = 1:size(WCenPos,2)
        spaDist(round(WCenPos{iPatt}(:,2)),round(WCenPos{iPatt}(:,1))) = ...
            spaDist(round(WCenPos{iPatt}(:,2)),round(WCenPos{iPatt}(:,1)))+1;
    end
    for iPatt = 1:size(WCenNeg,2)
        spaDist(round(WCenNeg{iPatt}(:,2)),round(WCenNeg{iPatt}(:,1))) = ...
            spaDist(round(WCenNeg{iPatt}(:,2)),round(WCenNeg{iPatt}(:,1)))+1;
    end
    
    % long
    for iPatt = 1:size(WCenPos,2)
        if size(WCenPos{iPatt},1) > 10
        spaDistlong(round(WCenPos{iPatt}(:,2)),round(WCenPos{iPatt}(:,1))) = ...
            spaDistlong(round(WCenPos{iPatt}(:,2)),round(WCenPos{iPatt}(:,1)))+1;
    
        end
    end
    for iPatt = 1:size(WCenNeg,2)
        if size(WCenNeg{iPatt},1) > 10
        spaDistlong(round(WCenNeg{iPatt}(:,2)),round(WCenNeg{iPatt}(:,1))) = ...
            spaDistlong(round(WCenNeg{iPatt}(:,2)),round(WCenNeg{iPatt}(:,1)))+1;
    
        end
    end
  toc
  
% tic
% load([pwd,'/preprocess/Preprocess_eigen345_right/Preprocessed_sub',num2str(iSub),'.mat'])
% sigIn = squeeze(sigBPass) ;
% 
% flagFindPhase = 1 ;
% if flagFindPhase
%     tic
%     hilbertSig = [] ;
%     phaseSig = [] ;
%     
%     for iX = 1:size(sigIn,1)
%         for iY = 1:size(sigIn,2)
%             hilbertSig(iX,iY,:) = hilbert(sigIn(iX,iY,:)) ;
%             phaseSig(iX,iY,:) = angle(hilbertSig(iX,iY,:)) ;
%         end
%     end
%     vPhaX = zeros(size(phaseSig)) ;
%     vPhaY = zeros(size(phaseSig)) ;
%     for iTime2 = 1:size(phaseSig,3)
%         for iX = 1:size(phaseSig,1)
%             vPhaX(iX,2:end-1,iTime2) = anglesubtract(phaseSig(iX,3:end,iTime2),phaseSig(iX,1:end-2,iTime2))/2 ;
%         end
%         for iY = 1:size(phaseSig,2)
%             vPhaY(2:end-1,iY,iTime2) = anglesubtract(phaseSig(3:end,iY,iTime2),phaseSig(1:end-2,iY,iTime2))/2 ;
%         end
%     end
%     temp = vPhaX ;
%     vPhaseX = -vPhaX./sqrt(vPhaX.^2+vPhaY.^2) ;
%     vPhaseY = -vPhaY./sqrt(temp.^2+vPhaY.^2) ;
%     toc
% end
%     
% curlz = [] ;
% cav = [] ;
% for time = 1:size(sigIn,3)
%     temp1_vx = vPhaseX(:,:,time);
%     temp1_vy = vPhaseY(:,:,time);
%     [x,y] = meshgrid(1:size(temp1_vx,2),1:size(temp1_vx,1));
%     
%     [curlz(:,:,time) cav(:,:,time)]  = curl(x,y,temp1_vx,temp1_vy);
% end
% cav_filt_pos = curlz;
% cav_filt_pos(cav_filt_pos<paramsVor.curlThre) = 0;
% %         cav_filt(cav_filt>0.6) = 1;
% cav_filt_pos(isnan(cav_filt_pos)) = 0;
% 
% cav_filt_neg = curlz;
% cav_filt_neg(cav_filt_neg>-paramsVor.curlThre) = 0;
% %         cav_filt(cav_filt>0.6) = 1;
% cav_filt_neg(isnan(cav_filt_neg)) = 0;
% 
% spaVorDist = spaVorDist+sum(cav_filt_pos,3)-sum(cav_filt_neg,3) ;

    disp(['finishing ...',num2str(iSub)])
toc
end

%%
% save('vorSpaceDist.mat','spaDist','spaVorDist','spaFullVorDist')
figure;
% imagesc(spaDist)
c = parula(256) ; c(1:10,:) = ones(10,3) ;
imagesc(spaDistlong)
colormap(c)
% imagesc(spaVorDist)
% imagesc(spaFullVorDist)

set(gca,'yDir','normal')
% % 14 sec
% label = load([pwd,'/Data/fMRI/network_array.csv']) ;  % 7 net for hierarchy
% 
% parLabel = ft_read_cifti([pwd,'/Data/fMRI2/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil_Colors.32k_fs_LR.dlabel.nii']) ;
% parValid2 = parLabel ;
% parValid2.x1(isnan(parLabel.x1),:) = [] ;
% posValid = posLC ;
% posValid.vertices(isnan(parLabel.x1) ,:) = [] ;
% 
% parValid = size(parValid2.x1) ;
% for iPar = min(parValid2.x1):max(parValid2.x1)
%     idxPar = find(parValid2.x1==iPar);
%     parValid(idxPar) = label(iPar);
% end

% load HCP 22 sections
addpath(genpath([pwd,'/ToolOthers/fMRI/cifti-matlab-master']))

folderName = [pwd,'/Data/fMRI2/100610_3T_structural/MNINonLinear/fsaverage_LR32k/'] ;
% fileName = '100610.L.flat.32k_fs_LR.surf.gii' ;   % left
fileName = '100610.R.flat.32k_fs_LR.surf.gii' ;   % right
posFile = [folderName,fileName] ;
posLC = gifti(posFile) ;

% 22 sec
parLabel2 = ft_read_cifti([pwd,'/Data/HCP/ResultsRegions_ROI.dlabel.nii']) ;
sec = 32493:32492*2 ; %right    (1:32492) %left
parValid = parLabel2.x1(sec) ;
parValid(isnan(parLabel2.x1(sec)),:) = [] ;
posValid = posLC ;
posValid.vertices(isnan(parLabel2.x1(sec)) ,:) = [] ;
for iPar = (1:22)
    y = [] ;
    x = [] ;
    [k_full,y_full,x_full] = genBoundary22(iPar,parValid,posValid,params) ;
    % idxA4 = find(parValid==iPar);
    for iSec = 1:length(y_full)
        y = [y;y_full{iSec}(k_full{iSec}) ] ;
        x = [x;x_full{iSec}(k_full{iSec}) ] ;
        hold on
        plot(y_full{iSec}(k_full{iSec}),x_full{iSec}(k_full{iSec}),'.','color',[0.8, 0, 0],'markerSize',2)
    end
end
 load([pwd,'/Results_data/fMRI_data/maskInfo_right.mat'])
    b_plot = [(b(:,1)-min(params.xCord)+1)/2+0.5, (b(:,2)-min(params.yCord)+1)/2+0.5 ] ;
    lh = plot(b_plot(:,1),b_plot(:,2),'r.','linewidth',0.5,'markerSize',2) ; %0.5)
    lh.Color = [lh.Color 0.4] ;
    axis off
    exportgraphics(gcf,['colorBarDist.pdf'],'ContentType','vector')

%% vortex trajectories
iSub = 4 ;
load([pwd,'/preprocess/VortexStat_eigen345/Vortex_sub',num2str(iSub),'.mat'])

DurPos = [] ;
for iPatt = 1:size(WCenPos,2)
    DurPos(iPatt) = size(WCenPos{iPatt},1) ;
end
DurNeg = [] ;
for iPatt = 1:size(WCenNeg,2)
    DurNeg(iPatt) = size(WCenNeg{iPatt},1) ;
end

DurMean = [] ;
for iTime = 1:1200
    idxPos = find(iTime>=pattStPos &iTime<=pattStPos+...
    DurPos-1) ; 
    idxNeg = find(iTime>=pattStNeg &iTime<=pattStNeg+...
    DurNeg-1) ; 
    DurMean(iTime) = mean([DurPos(idxPos),DurNeg(idxNeg) ]) ;
    
end
[~,sortIdx] = sort(DurMean) ;

[~,sortIdx2] = sort([DurNeg,DurPos]) ;

%% 
iTime = sortIdx(end-60) ;   % 150,25,
pattStFull = [pattStNeg,pattStPos] ; iTime = pattStFull(sortIdx2(end-10)) ;
figure;
n=30;
% cmap = [linspace(0,0.9,n)', linspace(.9447,.447,n)', linspace(.9741,.741,n)'];
cmap = [linspace(0.99,0.99,n)', linspace(.0,1.0,n)', linspace(.0,.5,n)'];
paramTraj.c = cmap ;
% paramTraj.c = [1,0,0] ;
paramTraj.linWid = 4 ;
% cmap = [ linspace(.9741,.741,n)', linspace(.9447,.447,n)',linspace(0,0.9,n)'];
cmap = [linspace(.0,.5,n)', linspace(.0,1.0,n)',linspace(0.99,0.99,n)'];
paramTraj2.c = cmap ;
% paramTraj2.c = [0,0,1] ;
paramTraj2.linWid = 4 ;
visTraj3(WCenPos, pattStPos,DurPos,iTime, paramTraj)

visTraj3(WCenNeg, pattStNeg,DurNeg,iTime, paramTraj2)

% load HCP 22 sections
addpath(genpath([pwd,'/ToolOthers/fMRI/cifti-matlab-master']))

folderName = [pwd,'/Data/fMRI2/100610_3T_structural/MNINonLinear/fsaverage_LR32k/'] ;
fileName = '100610.L.flat.32k_fs_LR.surf.gii' ;
posFile = [folderName,fileName] ;
posLC = gifti(posFile) ;

% 22 sec
parLabel2 = ft_read_cifti([pwd,'/Data/HCP/ResultsRegions_ROI.dlabel.nii']) ;
parValid = parLabel2.x1(1:32492) ;
parValid(isnan(parLabel2.x1(1:32492)),:) = [] ;
posValid = posLC ;
posValid.vertices(isnan(parLabel2.x1(1:32492)) ,:) = [] ;
for iPar = 1:22%(1:22)
    y = [] ;
    x = [] ;
    [k_full,y_full,x_full] = genBoundary22(iPar,parValid,posValid,params) ;
    % idxA4 = find(parValid==iPar);
    for iSec = 1:length(y_full)
        y = [y;y_full{iSec}(k_full{iSec}) ] ;
        x = [x;x_full{iSec}(k_full{iSec}) ] ;
        plot(y_full{iSec}(k_full{iSec}),x_full{iSec}(k_full{iSec}),'color',[0.6, 0.6, 0.6])
    end
end
saveName = 'wholeBrain' ;
% exportgraphics(gcf,['vortices_',saveName,'.pdf'],'ContentType','vector')
% exportgraphics(gcf,['traj2LP','.pdf'],'ContentType','vector')

% figure;
% colormap(cmap)
% colorbar

%% fig. 8A. schematic vortex array (draw on illustrator)



%% speed
v = [] ;
for iSub = 1:180
    
for iBurst = 1:size(WCenNeg,2)
    for iStep = 1:size(WCenNeg{iBurst},1)
        v = [v;sqrt(diff(WCenNeg{iBurst}(:,1)).^2+ diff(WCenNeg{iBurst}(:,2)).^2)/0.72] ;
    end
end

end

%% supplementary movie
visX = 96:124 ; visY = 168:196 ;
vDownR = 1 ;
visX = 1:176 ; visY = 1:251 ;
vDownR = 3 ;

flagSaveVideo = 1 ;
if flagSaveVideo
    t =  datetime('now') ;
    dateStr = datestr(t,'mmmmdd_HH:MM') ;
    vidTitle = [pwd,'/Results/Project3/movies/'] ;
    saveFileName = [vidTitle,'fMRI_s1','_sup_movie1_3_',dateStr] ;
    vidObj = VideoWriter(saveFileName,'Motion JPEG AVI');
    v.Quality = 50 ;
    vidObj.FrameRate = 5 ;  % 5
    open(vidObj);
    fig = figure;
    set(gcf,'Position', [85 192 1064 420])
    % set(gcf,'Position',[358 266 560 420])
    % set(gcf,'Position',[39 251 1168 435]) % for 2 subplot
else
    figure;
    set(gcf,'Position', [85 192 1064 420])
end

timeCount = 0 ;
for iTime = 350:449
    timeCount = timeCount+1 ;
subplot(1,2,1)
phaseSig = phaseSig.*mask ;
imagesc(visY,visX,phaseSig(visX,visY,iTime))
set(gca,'yDir','normal')
caxis([-pi-0.1 pi])
% c = pmkmp_new ; 
c1 = jet(256) ;
c = c1(20:236,:) ;
% c = parula(256)
c(1,:) = 1 ;
colormap(c)
colorbar('southoutside')
hold on
% line([visY2(1),visY2(end)],[visX2(1),visX2(1)],'Color',[0.4,0.4,0.4],'LineStyle','--','LineWidth',1) ;
% line([visY2(1),visY2(end)],[visX2(end),visX2(end)],'Color',[0.4,0.4,0.4],'LineStyle','--','LineWidth',1) ;
% line([visY2(1),visY2(1)],[visX2(1),visX2(end)],'Color',[0.4,0.4,0.4],'LineStyle','--','LineWidth',1) ;
% line([visY2(end),visY2(end)],[visX2(1),visX2(end)],'Color',[0.4,0.4,0.4],'LineStyle','--','LineWidth',1) ;
hAx = gca;
set(hAx, 'box','on','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
title(['phase field of brain activity at timeframe ',num2str(timeCount)])
view( az , el )
axis off

subplot(1,2,2)
imagesc(visY,visX,cavPlot(visX,visY,iTime)/2.5)
set(gca,'yDir','normal')
caxis([-2 2])
m = 256 ;
if (mod(m,2) == 0)
    m1 = m*0.25;
    r = [zeros(m1,1); (0:m1-1)'/max(m1-1,1)] ;
    % r2 = [ones(m1,1); (m1-1:-1:0)'/max(m1-1,1)/2+0.5 ] ;
    r2 = [ones(m1,1); ones(m1,1)] ;
    r = [r; r2];
    g = [zeros(m1,1);(0:m1-1)'/max(m1-1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
    
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end
c2 = [r g b];
c2(1,:) = 1 ;
colormap(gca,c2)
colorbar('southoutside')
% caxis([-1 1])

hold on
% vDownR = 3 ;
quiver(visY(1:vDownR:end),visX(1:vDownR:end),vPhaseX(visX(1:vDownR:end),...
    visY(1:vDownR:end),iTime),vPhaseY(visX(1:vDownR:end),visY(1:vDownR:end),iTime)...
    ,0.5,'Color',[0.6 0.6 0.6], 'MaxHeadSize',1)

hold on
% line([visY2(1),visY2(end)],[visX2(1),visX2(1)],'Color',[0.6,0.6,0.6],'LineStyle','--','LineWidth',1) ;
% line([visY2(1),visY2(end)],[visX2(end),visX2(end)],'Color',[0.6,0.6,0.6],'LineStyle','--','LineWidth',1) ;
% line([visY2(1),visY2(1)],[visX2(1),visX2(end)],'Color',[0.6,0.6,0.6],'LineStyle','--','LineWidth',1) ;
% line([visY2(end),visY2(end)],[visX2(1),visX2(end)],'Color',[0.6,0.6,0.6],'LineStyle','--','LineWidth',1) ;

hAx2 = gca;
set(hAx2, 'box','on','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
title(['vorticities and detected vortices at timeframe ',num2str(timeCount)])

view( az , el )
axis off
toc

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