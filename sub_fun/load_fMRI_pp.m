function [sigBPass] = load_fMRI_pp(...
    dataDir,posFile,params)
% example pararms
% params.sigmScale = 2.^(6:-1:2) ; params.flagBP = 1;
% params.downSRate = 2 ;                        % downsample the re-interpolation
% params.xCord = -250:params.downSRate:250 ;    % coordinate re-interpolation
% params.yCord = -150:params.downSRate:200 ;
% params.fsTem = 1/0.72 ;                       % temporal sampling rate
% params.zscore = 0 ;                           % flag, 1 for zscore, 0 for no     
% params.flagSur =0  ; params.surMethodNum =3 ; 
% params.flagVisBpSig =0 ;

disp('start importing and pre-processing ...')
tic


downSRate = params.downSRate ;
xCord =params.xCord ;
yCord = params.yCord ;
flagSur = params.flagSur ; 
surMethodNum = params.surMethodNum ; 
flagVisBpSig = params.flagVisBpSig ;
flagBP = params.flagBP ;

%% Import cortex data
tic
% addpath(genpath([pwd,'/ToolOthers/fMRI/cifti-matlab-master']))
% load the surface data
data = ft_read_cifti(dataDir) ;
% load the position data
posLC = gifti(posFile) ;

fsTem = 1/0.72 ;
% this function extract the cortex part of data from the HCP data
[sigValid,posValid,nanChans] = preproc_fRMI(data,posLC,fsTem) ;
toc
clearvars data posLC posRC

% checkRegion(posValid)

%% new interpolation based on the whole left cortex
x = double(posValid.vertices(:,1));
y = double(posValid.vertices(:,2));
% k = convhull(x,y);
k = alphaShape(x,y,4) ;
[a, b] = k.boundaryFacets();
figure;
plot(x,y,'bo')
hold on;
plot(b(:,1),b(:,2),'r','linewidth',8)
clearvars x y k a

%
% downSRate = 2 ;  % 1 or 2 or 3
% xCord = -250:downSRate:250 ;
% yCord = -150:downSRate:200 ;

% find the mask by round original points (but sparse)
% temp = zeros(length(xCord),length(yCord)) ;
% indValid = sub2ind(size(temp),round(posValid.vertices(:,1))-min(xCord)+1,...
%     round(posValid.vertices(:,2))-min(yCord)+1) ;
% temp(indValid) = 1 ;
% imagesc(temp)

% find the mask by filling the boundary
bw = poly2mask(b(:,1)-min(xCord)+1,b(:,2)-min(yCord)+1,max(yCord)-min(yCord)+1,...
    max(xCord)-min(xCord)+1);
% figure;
% imshow(bw)
% hold on
% plot(b(:,1)-min(xCord)+1,b(:,2)-min(yCord)+1,'b','LineWidth',2)
% hold off
        
% disp(['new boundary after interpolation: x: ',num2str(min(xCord)), ', '...
%     num2str(max(xCord)),' y:' num2str(min(yCord)), ',',...
%     num2str(max(yCord)) ])
flagPlotSpe = 0 ;
flagCheckIntp = 0 ;

% sigOri_reInterp = spaceFreq_fMRI(posValid,sigValid,xCord,yCord,...
%     flagPlotSpe, flagCheckIntp) ;
mask = double(bw(1:downSRate:end,1:downSRate:end)) ;
mask(mask==0) = nan ;
data_reInterp = spaceFreq_fMRI(posValid,sigValid,xCord,yCord,...
    flagPlotSpe, flagCheckIntp,mask) ;
clearvars a bw flagPlotSpe flagCheckIntp


%% generate surrogate data
% flagSur = 0 ;
if flagSur
    disp('surrogating data ...')
    % surMethodNum = 4 ;
    dataOut = surrogate_fMRI(data_reInterp, [], surMethodNum) ;
    dataOut = dataOut.*mask ;

else
    disp('no surrogate ...')
    dataOut = data_reInterp ;
end

%% bandpass filtering (before re-interpolation)
close all
fLow = 0.01 ;
fHigh = 0.1 ;
dataOut_reshape = reshape(dataOut,size(dataOut,1)*size(dataOut,2),[]) ;
[bandpasSig,~, ~, ~] = bandpa_fMRI(dataOut_reshape,fsTem,fLow,fHigh) ;
if params.zscore == 1
    bandpasSig = zscore(bandpasSig,[],2) ;
end
bandpass_reInterp = reshape(bandpasSig,size(dataOut,1),size(dataOut,2),[]) ;

%% spatial filter 1 - Gaussian (after re-interpolation)
% sigmScale = [8,6,4,2,0.1] ;
% sigmScale = 2.^(6:-2:-2) ;
% sigmScale = 2.^(3:-1:-1) ;
% sigmScale = 2.^(6:-1:2)/downSRate ;
sigmScale = params.sigmScale/downSRate ;
if flagBP==1
    numScale = length(sigmScale)-1 ; 
else
    numScale = length(sigmScale) ;
end
sigLPass = [] ;
sigBPass = [] ;

%% bandpass filter
if flagBP==1
for iTime = 1:size(bandpass_reInterp,3)
%     sigLPass(1,:,:,iTime) = imgaussfilt(bandpass_reInterp(:,:,iTime),sigmScale(1)) ;
    filtSigma = sigmScale(1);   % 0.6
    filtWidth = ceil(3*filtSigma);
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    sigLPass(1,:,:,iTime) = nanconv(bandpass_reInterp(:,:,iTime),imageFilter,'edge', 'nanout');

    for iScale = 1:numScale     
        % without nan
        % sigLPass(iScale+1,:,:,iTime) = imgaussfilt(bandpass_reInterp(:,:,iTime),sigmScale(iScale+1)) ;
        % with nan
        filtSigma = sigmScale(iScale+1);   % 0.6
        filtWidth = ceil(3*filtSigma);
        imageFilter=fspecial('gaussian',filtWidth,filtSigma);
        sigLPass(iScale+1,:,:,iTime) = nanconv(bandpass_reInterp(:,:,iTime),imageFilter,'edge', 'nanout');
        sigBPass(iScale,:,:,iTime) = sigLPass(iScale+1,:,:,iTime) - sigLPass(iScale,:,:,iTime) ;
    end
end
disp('total time consumption is...')
toc

elseif flagBP==2  % low pass
    for iTime = 1:size(bandpass_reInterp,3)
%     sigLPass(1,:,:,iTime) = imgaussfilt(bandpass_reInterp(:,:,iTime),sigmScale(1)) ;
    for iScale = 1:numScale     
        % without nan
        % sigLPass(iScale+1,:,:,iTime) = imgaussfilt(bandpass_reInterp(:,:,iTime),sigmScale(iScale+1)) ;
        % with nan
        filtSigma = sigmScale(iScale);   % 0.6
        filtWidth = ceil(3*filtSigma);
        imageFilter=fspecial('gaussian',filtWidth,filtSigma);
        sigBPass(iScale,:,:,iTime) = nanconv(bandpass_reInterp(:,:,iTime),imageFilter,'edge', 'nanout');
    end
    end
    
    
%% low pass filter
else              % LoG
for iTime = 1:size(bandpass_reInterp,3)
%     sigLPass(1,:,:,iTime) = imgaussfilt(bandpass_reInterp(:,:,iTime),sigmScale(1)) ;
    for iScale = 1:numScale     
        % without nan
        % sigLPass(iScale+1,:,:,iTime) = imgaussfilt(bandpass_reInterp(:,:,iTime),sigmScale(iScale+1)) ;
        % with nan
        filtSigma = sigmScale(iScale);   % 0.6
        filtWidth = ceil(3*filtSigma);
        imageFilter=fspecial('log',filtWidth,filtSigma);
        sigBPass(iScale,:,:,iTime) = nanconv(bandpass_reInterp(:,:,iTime),imageFilter,'edge', 'nanout');
    end
end
disp('total time consumption is...')
% toc
end
%%
if flagVisBpSig
    iTime = 10 ;
    figure;
    imagesc(squeeze(bandpass_reInterp(:,:,iTime)))
    for iScale = 1:numScale
    % figure;
    % imagesc(squeeze(sigLPass(iScale,:,:,iTime)))
    figure;
    imagesc(squeeze(sigBPass(iScale,:,:,iTime)))
    end
end