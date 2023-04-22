function [dataOut, b, bandpass_reInterp,data_reInterp,posValid,nanChans,mask] = load_fMRI(...
    dataDir,posFile,flagSur,surMethodNum,params,flagVisBpSig,flagRest,flagSmooth,hemisphere,flagTask)

disp('start importing and pre-processing ...')

downSRate = params.downSRate ;
xCord =params.xCord ;
yCord = params.yCord ;

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
% clearvars data posLC posRC

% checkRegion(posValid)

%% new interpolation based on the whole left cortex
x = double(posValid.vertices(:,1));
y = double(posValid.vertices(:,2));
% k = convhull(x,y);
k = alphaShape(x,y,4) ;
[a, b] = k.boundaryFacets();

bw = poly2mask(b(:,1)-min(xCord)+1,b(:,2)-min(yCord)+1,max(yCord)-min(yCord)+1,...
    max(xCord)-min(xCord)+1);

flagPlotSpe = 0 ;
flagCheckIntp = 0 ;

mask = double(bw(1:downSRate:end,1:downSRate:end)) ;
mask(mask==0) = nan ;
data_reInterp = spaceFreq_fMRI(posValid,sigValid,xCord,yCord,...
    flagPlotSpe, flagCheckIntp,mask) ;
dataOut = data_reInterp;
clearvars a bw flagPlotSpe flagCheckIntp
%% temporal bandpass filtering 
if flagSmooth == 1 || flagSmooth == 2
    close all
    dataIn = dataOut;
    fLow = 0.01 ; % lower limit of the bandpass filter
    fHigh = 0.1 ; % upper limit of the bandpass filter
    dataOut_reshape = reshape(dataIn,size(dataIn,1)*size(dataIn,2),[]) ;
    [bandpasSig,~, ~, ~] = bandpa_fMRI(dataOut_reshape,fsTem,fLow,fHigh) ;
%     bandpasSig = zscore(bandpasSig,[],2) ; % ***************
    bandpass_reInterp = reshape(bandpasSig,size(dataIn,1),size(dataIn,2),[]) ;
    dataOut = bandpass_reInterp;
end
%% generate surrogate data
% if flagSur == 1
%     % 3D fourier transform-based random phase shuffling
%     if flagRest == 0 % task data
%         dataIn = dataOut(1:175,1:251,1:315); % ensure odd number
%     elseif flagRest == 1 % rest data
%         dataIn = dataOut(1:175,1:251,1:1199); % rest data
%     end
%         surData = zeros(size(dataIn)) ; 
%         dataIn(isnan(dataIn)) = 0 ;  
%         % find the middle point
%         midPoint = floor(size(dataIn)/2)+1 ;        
%         indMid = sub2ind(size(dataIn),midPoint(1),midPoint(2),midPoint(3)) ;
%  
%             randSeq = randperm(indMid - 1) ; % randomization
%             freqData = fftshift(fftn(dataIn(:,:,:))) ; % 3D fourier transform ??????
%             phaseData = angle(freqData) ;
%             absData = abs(freqData) ;
%             % randomize all the phase before the mid point (guarantee symmetry/mirror)
%             phaseData(1:indMid-1) = rand(indMid-1,1).*2*pi-pi; % random phase values -pi~pi
%             for iInd = 1:indMid - 2
%                 phaseData(indMid+iInd) = -phaseData(indMid-iInd) ; 
%             end
%             freqDataNew = ifftshift(absData.*exp(1i*phaseData)) ;
%             surData(:,:,:) = real(ifftn(freqDataNew)) ; % 3D inverse fourier tranform
% %         end
%          dataOut = surData;
%          
% 
% else
%     disp('no surrogate ...')
% end

if flagSur == 1
   if flagRest == 0 && flagTask ~= 3 % task data
        dataIn = dataOut(1:175,1:251,1:315); % ensure odd number
    elseif flagRest == 1 % rest data
        dataIn = dataOut(1:175,1:251,1:1199); 
   end
   if flagRest == 0 && flagTask == 3
      dataIn = dataOut(1:175,1:251,1:405); % ensure odd number
   end

    surData = zeros(size(dataIn)) ; 
    dataIn(isnan(dataIn)) = 0 ;  
    midPoint = floor(size(dataIn)/2)+1 ;        
    indMid = sub2ind(size(dataIn),midPoint(1),midPoint(2),midPoint(3)) ;
    randSeq = randperm(indMid - 1) ; 
    freqData = fftshift(fftn(dataIn)) ; % 3D fourier transform 
    phaseData = angle(freqData) ;
    absData = abs(freqData) ;
    phaseDataRand = phaseData;
    phaseDataRand(1:indMid-1) = rand(indMid-1,1).*2*pi-pi; % random phase values -pi~pi
    for iInd = 1:indMid - 2
        phaseDataRand(indMid+iInd) = -phaseDataRand(indMid-iInd) ; 
    end
    freqDataNew = ifftshift(absData.*exp(1i.*phaseData + 1i*phaseDataRand)) ; % random phase shuflling
    dataOut = real(ifftn(freqDataNew)) ; % 3D inverse fourier tranform
    dataOut_unsmooth = dataOut;
end  

    %% spatial bandpass filtering - Gaussian 
    if flagSur == 0 && flagSmooth == 2 
        dataIn = dataOut(1:175,:,:);
        sigmScale = params.sigmScale/downSRate ;
        sigLPass = [] ;
        sigBPass = [] ;
        % spatial bandpass filter
        numScale = length(sigmScale)-1 ;
        for iTime = 1:size(dataIn,3)
            filtSigma = sigmScale(1);   % 0.6
            filtWidth = ceil(3*filtSigma);
            imageFilter=fspecial('gaussian',filtWidth,filtSigma);
            sigLPass(1,:,:,iTime) = nanconv(dataIn(:,:,iTime),imageFilter,'edge', 'nanout');
            for iScale = 1:numScale     
                filtSigma = sigmScale(iScale+1);  
                filtWidth = ceil(3*filtSigma);
                imageFilter=fspecial('gaussian',filtWidth,filtSigma);
                sigLPass(iScale+1,:,:,iTime) = nanconv(dataIn(:,:,iTime),imageFilter,'edge', 'nanout');
                sigBPass(iScale,:,:,iTime) = sigLPass(iScale+1,:,:,iTime) - sigLPass(iScale,:,:,iTime) ;
            end
        end
        
        dataOut = permute(sigBPass,[2,3,4,1]);
        
    elseif  flagSur == 1 && flagSmooth == 1 
        dataIn = dataOut(1:175,:,:);
        sigmScale = params.sigmScale/downSRate ;
        sigLPass = [] ;
        sigBPass = [] ;
        % spatial bandpass filter
        numScale = length(sigmScale)-1 ;
        for iTime = 1:size(dataIn,3)
            filtSigma = sigmScale(1);   % 0.6
            filtWidth = ceil(3*filtSigma);
            imageFilter=fspecial('gaussian',filtWidth,filtSigma);
            sigLPass(1,:,:,iTime) = nanconv(dataIn(:,:,iTime),imageFilter,'edge', 'nanout');
            for iScale = 1:numScale     
                filtSigma = sigmScale(iScale+1);  
                filtWidth = ceil(3*filtSigma);
                imageFilter=fspecial('gaussian',filtWidth,filtSigma);
                sigLPass(iScale+1,:,:,iTime) = nanconv(dataIn(:,:,iTime),imageFilter,'edge', 'nanout');
                sigBPass(iScale,:,:,iTime) = sigLPass(iScale+1,:,:,iTime) - sigLPass(iScale,:,:,iTime) ;
            end
        end
        
        dataOut = permute(sigBPass,[2,3,4,1]);
    end
    
    
    %% ensure the edge of the output data is compatible with the flattenedd cortex
    if hemisphere == 1
        load('parcellation_template.mat')
    elseif hemisphere == 2
        load('parcellation_template22_RightBrain_subject1-100.mat')
        parcellation_template = parcellation_template22_RightBrain_100sub(:,:,1);
    end
    dataOut = dataOut(1:175,:,:).*parcellation_template(1:175,:)./parcellation_template(1:175,:);
    if flagSur == 1
       dataOut_unsmooth = dataOut_unsmooth.*parcellation_template(1:175,:)./parcellation_template(1:175,:);
       dataOut(:,:,:,1) = dataOut;
       dataOut(:,:,:,2) = dataOut_unsmooth; % record both smooth and unsmoothed surrogate data in the same output
    end
end
