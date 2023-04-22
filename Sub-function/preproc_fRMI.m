function [sigValid,posValid,nanChans] = preproc_fRMI(data,pos,fs) ;
% load in the signals from HCP and extract the valid surface cortex data
% input: data from nii file;
%        pos  from surf.nii file;
% output: sigValid and posValid, valid cortex data

disp('start preprocessing...');
surfData = data.dtseries ;
surfPos = data.pos ;

% cortex position
% for i=1:21
%     data.brainstructurelabel{i} 
% end
% cortex is the region we are interested but the pos is from surface.gii
% file

% find the cortex
cortexLIdx = find(data.brainstructure == 1) ;
startL = min(cortexLIdx) ;
endL = max(cortexLIdx) ;
cortexRIdx = find(data.brainstructure == 2) ;
startR = min(cortexRIdx) ;
endR = max(cortexRIdx) ;

%% Get the valid region (cortex)
%patch('Faces',posLC.faces,'Vertices',posLC.vertices)
%axis equal
rangeC = startL:endL ;  % left brain
% rangeC = startR:endR ;  % right brain

tic
% check the valid region in time
% figure; 
% plot((1:size(surfData,2))/fs,surfData(startL,:))
% hold on
% plot((1:size(surfData,2))/fs,surfData(startL+1000,:))
% hold on
% plot((1:size(surfData,2))/fs,surfData(startL+10000,:))
% title('check three lines of data to confirm the valid region')

% get rid of the first and the last 10s
% timeRange = fix(10*fs): size(surfData,2)-fix(10*fs) ;
timeRange = 1: size(surfData,2) ;
sigRange = surfData(rangeC,:) ;
sigOri = sigRange(:,timeRange) ;
% sigOri = zscore(surfData(rangeC,timeRange),[],2) ;
nanChans = any(isnan(sigOri(:,:)),2);

sigValid = sigOri;
sigValid(nanChans,:) = [] ;
posValid = pos ;
posValid.vertices(nanChans,:) = [] ;

% figure;
% plot(posValid.vertices(1:10:end,1),posValid.vertices(1:10:end,2),'o')