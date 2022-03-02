%% spatial power spectrum

params.downSRate = 2 ;                        % downsample the re-interpolation
params.xCord = -250:params.downSRate:250 ;    % coordinate re-interpolation
params.yCord = -150:params.downSRate:200 ;
params.fsTem = 1/0.72 ;                       % temporal sampling rate
params.zscore = 0 ;                           % flag, 1 for zscore, 0 for no     
params.flagSur =0  ; params.surMethodNum =3 ; 
params.flagVisBpSig =0 ; 

iSub = 1 ;

name = dir([pwd,'/Data/fMRI_converted_language']) ;
    dataDir = [name(iSub+5).folder,'/',name(iSub+5).name,'/tfMRI_LANGUAGE_LR_Atlas.dtseries.nii'] ;
name = dir([pwd,'/Data/fMRI_Converted']) ;
    posFile = [name(iSub+2).folder,'/',name(iSub+2).name,'/L.flat.32k_fs_LR.surf.gii'] ;

downSRate = params.downSRate ;
xCord =params.xCord ;
yCord = params.yCord ;
flagSur = params.flagSur ; 
surMethodNum = params.surMethodNum ; 
flagVisBpSig = params.flagVisBpSig ;
flagBP = params.flagBP ;

% Import cortex data
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

% new interpolation based on the whole left cortex
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

% find the mask by filling the boundary
bw = poly2mask(b(:,1)-min(xCord)+1,b(:,2)-min(yCord)+1,max(yCord)-min(yCord)+1,...
    max(xCord)-min(xCord)+1);

flagPlotSpe = 0 ;
flagCheckIntp = 0 ;

% sigOri_reInterp = spaceFreq_fMRI(posValid,sigValid,xCord,yCord,...
%     flagPlotSpe, flagCheckIntp) ;
mask = double(bw(1:downSRate:end,1:downSRate:end)) ;
mask(mask==0) = nan ;
data_reInterp = spaceFreq_fMRI(posValid,sigValid,xCord,yCord,...
    flagPlotSpe, flagCheckIntp,mask) ;
clearvars a bw flagPlotSpe flagCheckIntp

%%

params.downSRate = 2 ;                        % downsample the re-interpolation
params.xCord = -250:params.downSRate:250 ;    % coordinate re-interpolation
params.yCord = -150:params.downSRate:200 ;
params.fsTem = 1/0.72 ;                       % temporal sampling rate
params.zscore = 0 ;                           % flag, 1 for zscore, 0 for no     
params.flagSur =0  ; params.surMethodNum =3 ; 
params.flagVisBpSig =0 ; 

iSub = 1 ;

name = dir([pwd,'/Data/fMRI_Converted']) ;
    dataDir = [name(iSub+2).folder,'/',name(iSub+2).name,'/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii'] ;
    posFile = [name(iSub+2).folder,'/',name(iSub+2).name,'/L.flat.32k_fs_LR.surf.gii'] ;

downSRate = params.downSRate ;
xCord =params.xCord ;
yCord = params.yCord ;
flagSur = params.flagSur ; 
surMethodNum = params.surMethodNum ; 
flagVisBpSig = params.flagVisBpSig ;
flagBP = params.flagBP ;

% Import cortex data
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

% new interpolation based on the whole left cortex
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

% find the mask by filling the boundary
bw = poly2mask(b(:,1)-min(xCord)+1,b(:,2)-min(yCord)+1,max(yCord)-min(yCord)+1,...
    max(xCord)-min(xCord)+1);

flagPlotSpe = 0 ;
flagCheckIntp = 0 ;

% sigOri_reInterp = spaceFreq_fMRI(posValid,sigValid,xCord,yCord,...
%     flagPlotSpe, flagCheckIntp) ;
mask = double(bw(1:downSRate:end,1:downSRate:end)) ;
mask(mask==0) = nan ;
data_reInterp2 = spaceFreq_fMRI(posValid,sigValid,xCord,yCord,...
    flagPlotSpe, flagCheckIntp,mask) ;
clearvars a bw flagPlotSpe flagCheckIntp

%%
close all
nfft = 1000 ;
dataFFTFull = zeros(nfft/2+1,1) ;
for iTime = 45%1:100 
    dataTest = reshape(squeeze(data_reInterp(:,:,iTime)),1,[]) ;
    dataTest(isnan(dataTest)) = [] ;
    dataTest = dataTest-mean(dataTest) ;
    % dataFFT = fft(dataTest,nfft) ;
    dataFFT = pmtm(dataTest,[],nfft) ;
    dataFFTFull = dataFFTFull+dataFFT(1:nfft/2+1) ;
end

dataFFTFull2 = zeros(nfft/2+1,1) ;
for iTime = 45%1:100 
    dataTest = reshape(squeeze(data_reInterp2(:,:,iTime)),1,[]) ;
    dataTest(isnan(dataTest)) = [] ;
    dataTest = dataTest-mean(dataTest) ;
    % dataFFT = fft(dataTest,nfft) ;
    dataFFT = pmtm(dataTest,[],nfft) ;
    dataFFTFull2 = dataFFTFull2+dataFFT(1:nfft/2+1) ;
end

figure; loglog(linspace(0,1,nfft/2+1),abs(dataFFTFull),'b')
hold on; loglog(linspace(0,1,nfft/2+1),abs(dataFFTFull2),'r')