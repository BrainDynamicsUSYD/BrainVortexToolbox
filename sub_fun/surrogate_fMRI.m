function surData = surrogate_fMRI(dataIn, posData, methodNum) 
% posData stores different data point on different rows and columns denotes
% different coordinates.
%%
switch methodNum
    case 1
        % keep spatial autocorrelation (SA) for arbitary co-ordinate
        % step 1, randomly shuffle the data
        posDataSur = posData(randperm(size(posData,1)),:) ;
        % step 2, smoothing and re-scaling
        for k = 10:10:100
            for iPoint = 1:size(posDataSur,1)
                [idx,d] = knnsearch(posDataSur,posDataSur(iPoint,:),'K',k) ;
                surData(iPoint) = kernelSm(d',dataIn(idx)) ;
            end
            disp(['finishing ', num2str(k),' %'])
        end
    case 2
        % keep SA by spatial FFT
        
    case 3
        % random shuffle
        surData = dataIn ;
        temp2 = surData (~isnan(surData)) ;   % get the valid data
        temp2 = temp2(randperm(numel(temp2))) ;          % please check the sequence of temp2 is changed
        surData(~isnan(surData)) = temp2 ;             
    
    case 4
        % fourier transport (regular coordinate)
        % freqData = zeros(size(dataIn)) ;
        % phaseData = zeros(size(dataIn)) ;
        % ampData = zeros(size(dataIn)) ;
        surData = zeros(size(dataIn)) ;
        dataIn(isnan(dataIn)) = 0 ;
        midPoint = floor(size(dataIn)/2)+1 ;
        for iTime = 1: size(dataIn,3)           
            freqData = fftshift(fft2(dataIn(:,:,iTime))) ;
            phaseData = angle(freqData) ;
            absData = abs(freqData) ;
            % randomize all the phase before the mid point (guarantee symmetry)
            indMid = sub2ind(size(phaseData),midPoint(1),midPoint(2)) ;
            phaseData(1:indMid-1) = phaseData(randperm(indMid-1)) ;
            for iInd = 1:indMid - 2
                phaseData(indMid+iInd) = -phaseData(indMid-iInd) ;
            end
            freqDataNew = ifftshift(absData.*exp(i*phaseData)) ;
            surData(:,:,iTime) = real(ifft2(freqDataNew)) ;
        end
    case 5
        % fourier transport (regular coordinate)
        % freqData = zeros(size(dataIn)) ;
        % phaseData = zeros(size(dataIn)) ;
        % ampData = zeros(size(dataIn)) ;
        surData = zeros(size(dataIn)) ;
        dataIn(isnan(dataIn)) = 0 ;
        midPoint = floor(size(dataIn)/2)+1 ;
        indMid = sub2ind(size(dataIn(:,:,1)),midPoint(1),midPoint(2)) ;
        randSeq = randperm(indMid - 1) ;
        for iTime = 1: size(dataIn,3)           
            freqData = fftshift(fft2(dataIn(:,:,iTime))) ;
            phaseData = angle(freqData) ;
            absData = abs(freqData) ;
            % randomize all the phase before the mid point (guarantee symmetry)
            
            phaseData(1:indMid-1) = phaseData(randSeq) ;
            for iInd = 1:indMid - 2
                phaseData(indMid+iInd) = -phaseData(indMid-iInd) ;
            end
            freqDataNew = ifftshift(absData.*exp(i*phaseData)) ;
            surData(:,:,iTime) = real(ifft2(freqDataNew)) ;
        end
    case 6
        
        surData = zeros(size(dataIn)) ;
        
        % temporal surrogate
        sigOriTemp = reshape(dataIn,size(dataIn,1)*size(dataIn,2),[]) ;
        surSigTemp = nan(size(sigOriTemp)) ;
        empChannels = find(isnan(sigOriTemp)) ;
        
        randNumHalf =  2*pi*rand(size(dataIn,3)/2-1,1) ;
        randNum = [0;randNumHalf;0;-flip(randNumHalf) ]' ;
        
        for iChannel = setdiff(1:size(sigOriTemp,1),empChannels)
            freqSig = fft(sigOriTemp(iChannel,:)) ;
            absFreq = abs(freqSig) ;
            phaseFreq = angle(freqSig) ;
            reconSig = ifft(absFreq.*exp(1i*phaseFreq)) ;
        
            surSigTemp(iChannel,:) = ifft(absFreq.*exp(1i*(phaseFreq+randNum))) ;
        end
        dataIn = reshape(real(surSigTemp),size(dataIn,1),size(dataIn,2),[]) ;
    
        % spatial surrogate
        surData = zeros(size(dataIn)) ;
%         dataIn(isnan(dataIn)) = 0 ;
%         midPoint = floor(size(dataIn)/2)+1 ;
%         for iTime = 1: size(dataIn,3)           
%             freqData = fftshift(fft2(dataIn(:,:,iTime))) ;
%             phaseData = angle(freqData) ;
%             absData = abs(freqData) ;
%             % randomize all the phase before the mid point (guarantee symmetry)
%             indMid = sub2ind(size(phaseData),midPoint(1),midPoint(2)) ;
%             phaseData(1:indMid-1) = phaseData(randperm(indMid-1)) ;
%             for iInd = 1:indMid - 2
%                 phaseData(indMid+iInd) = -phaseData(indMid-iInd) ;
%             end
%             freqDataNew = ifftshift(absData.*exp(i*phaseData)) ;
%             surData(:,:,iTime) = real(ifft2(freqDataNew)) ;
%         end

        dataIn(isnan(dataIn)) = 0 ;
        midPoint = floor(size(dataIn)/2) ;
        indMid = sub2ind(size(dataIn(:,:,1)),midPoint(1),midPoint(2)) ;

        for iTime = 1: size(dataIn,3)   
            randSeq = randperm(indMid - 1) ;
            freqData = fftshift(fft2(dataIn(:,:,iTime))) ;
            phaseData = angle(freqData) ;
            absData = abs(freqData) ;
            % randomize all the phase before the mid point (guarantee symmetry)
            
            phaseData(1:indMid-1) = phaseData(randSeq) ;
            for iInd = 1:indMid - 2
                phaseData(indMid+iInd) = -phaseData(indMid-iInd) ;
            end
            freqDataNew = ifftshift(absData.*exp(i*phaseData)) ;
            surData(:,:,iTime) = real(ifft2(freqDataNew)) ;
        end
end

end
%% smooth using the kernel function
function dataOut = kernelSm(distData,dataIn)
    dataOut = sum( exp(-distData).*dataIn )/sum( exp(-distData)) ;    
end

%%
function findVariogram(dataIn,posData)
    % find all distances
    
    % 25 uniformly spaced distance intervals
    
    % Smoothing variogram

end