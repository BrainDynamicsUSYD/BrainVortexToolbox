function [WCentroids,absoluteTime,instantTotalPower,pattSize,patternIdx]= pattDetection_v4(sigIn, sigBinary, params, flagSaveData,saveFileName)
% function for pattern detection
% in this version, 3D detection, max size of pattern are considered
%% spatial-temporal pattern detection based on continuity
% CC = bwconncomp(sigBinary,6) ;               % continued points in space-time
CC = bwconncomp(sigBinary) ;               % continued points in space-time

B1 = regionprops(CC,'BoundingBox');        % bounding box of patterns
boundPatt = cat(1, B1.BoundingBox);
areaPatt = regionprops(CC,'Area') ;            % for computing total scale of patterns


%% initialization
countPatt = 1;                             % for counting patterns
Duration = [] ;                            % pattern duration
pattSize = [] ;                            % pattern size
Centroids = [] ;                           % geometry centre of patterns
WCentroids = [] ;                          % weighted centre of patterns
distCent = [] ;                            % distance between patterns
centInterval = [] ;                        % intervel between patterns
instantScale = [] ;                        % instantaneous size in patterns
instantPeakAmp = [] ;                      % instantaneous peak amplitudue in patterns
instantTotalPower = [] ;                   % instantaneous total power in patterns
rangeFrame = [] ;                          % the start and end time frames of patterns
width = [] ;                               % instantaneous size of patterns


% params.minpatternTime = 30 ;             % for 1 Gamma cycle (30Hz), 3 cycles (80)
sigPlot = sigBinary.*sigIn ;               % for calculating weighted centres
fullBinary = zeros(size(sigPlot)) ;        % store the final pattern index

%% further pattern detection
for iPatt = 1: size(CC.PixelIdxList,2)
    currentIdx = CC.PixelIdxList{iPatt} ;        % extract index for patterns
    DurationAll(iPatt) = boundPatt(iPatt,6) ; 

    pattTimeStart = boundPatt(iPatt,3)+0.5 ;     % note that the first frame
                                                 % starts with 0.5
    pattTimeEnd = pattTimeStart + boundPatt(iPatt,6) -1 ;

    
    % temporal threshold
    if DurationAll(iPatt) < params.minPattTime
        continue
    end
    % spatial threshold
    if boundPatt(iPatt,4)< params.minPattSize || boundPatt(iPatt,5)<params.minPattSize
        continue
    end
    % Amp = [Amp; sigPlot(currentIdx)];
    
    % pattern properties to be stored
    Duration(countPatt) = pattTimeEnd-pattTimeStart+1 ;    % duration
    pattSize(countPatt) = areaPatt(iPatt).Area ;  % for total scale
    
    pattIdxTemp = zeros(size(sigPlot)) ;
    pattIdxTemp(currentIdx) = 1 ;
    currentPatt = sigPlot.*pattIdxTemp ;
    % sumAmp(count) = sum(currentPatt(:)) ;     % sum of amplitude
    % peakAmp(count) = max(currentPatt(:)) ;    % peak amplitude
    
    % loop through each time frame to study instaneous properties within
    % patterns
     absoluteTime{countPatt} = [pattTimeStart:pattTimeEnd]; % find time points for each pattern
     
    timeCount = 1 ;
    for iTime = pattTimeStart:pattTimeEnd
        % grab the current patterns
        instantBinary = pattIdxTemp(:,:,iTime) ;
        
        fullBinary(:,:,iTime) = fullBinary(:,:,iTime) + instantBinary ;
        instantPattern{countPatt}(:,:,timeCount) = instantBinary.*sigPlot(:,:,iTime) ;
        
        instantScale{countPatt}(timeCount,:) = sum(instantBinary(:)) ;  % instant scale
        
        tempPeak = max(max(currentPatt(:,:,iTime))) ;
        instantPeakAmp{countPatt}(timeCount,:) = tempPeak ;
        
        instantTotalPower{countPatt}(timeCount,:) = sum(sum(abs(currentPatt(:,:,iTime)) )) ;
        
       % find avg abs phase of electrodes in the pattern
        instantMeanPhase{countPatt}(timeCount,:) = sum(sum(abs(currentPatt(:,:,iTime))))./instantScale{countPatt}(timeCount,:);
       
        
        S = regionprops(instantBinary,instantPattern{countPatt}(:,:,timeCount),{'Centroid','WeightedCentroid'} );
        % calculate
        B2 = regionprops(instantBinary,'BoundingBox');
        
        boundPattInst = cat(1, B2.BoundingBox);
        width{countPatt}(timeCount,:) = (boundPattInst(3)+boundPattInst(4))/2 ;
        
        
        Centroids{countPatt}(timeCount,:) = cat(1, S.Centroid);
        WCentroids{countPatt}(timeCount,:) = cat(1, S.WeightedCentroid) ;
        timeCount = timeCount + 1;
    end
    rangeFrame(countPatt,:) = [pattTimeStart,pattTimeEnd] ;
    %
    if countPatt>1
        firstCentroidsLoc = squeeze(WCentroids{countPatt}(1,:)) ;
        % calculate the distance of centroids of two patterns
        distCent(countPatt) = sqrt(sum((firstCentroidsLoc - lastCentroidsLoc).^2)) ;
        
        firstCentroidsTime = pattTimeStart ;
        % calculate the time interval between patterns
        centInterval(countPatt) = firstCentroidsTime - lastCentroidsTime ;
    end
    lastCentroidsLoc = squeeze(WCentroids{countPatt}(end,:)) ;
    lastCentroidsTime = pattTimeEnd ;
    countPatt = countPatt+1 ;
end
CC_patterns = bwconncomp(fullBinary) ;       % continued points in space-time
patternIdx = CC_patterns.PixelIdxList ;      % save index to reduce size

%% saving data
    
%if flagSaveData
    
    %saveFileName = [saveFileName,'_v4.mat'] ;
    %save(saveFileName, 'Centroids','WCentroids', 'rangeFrame',...
    %    'pattSize','instantScale','Duration','DurationAll','distCent',...
    %    'centInterval','instantPeakAmp','instantTotalPower','width','patternIdx','instantPattern','absoluteTime','instantMeanPhase') ;
    

%end