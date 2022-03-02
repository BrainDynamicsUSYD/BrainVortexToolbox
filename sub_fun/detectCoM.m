function [durationValid,WCentroids,Centroids,pattStValid,areaPatt]= detectCoM(sigBinary,sigIn,paramCoM)
thresholdTime = paramCoM.thresholdTime ;
thresholdSize = paramCoM.thresholdSize ;

        % find the pattern
    CC = bwconncomp(sigBinary) ;               % continued points in space-time
    B1 = regionprops(CC,'BoundingBox');        % bounding box of patterns
    boundPatt = cat(1, B1.BoundingBox);
    area = regionprops(CC,'Area') ;
    areaPatt = cat(1, area.Area) ;
    
    
    tic
      % loop through each patterns to find instantaneous centre
    patCount = 1 ;
    tic
    for iPatt = 1: size(CC.PixelIdxList,2)
        if mod(iPatt,floor(size(CC.PixelIdxList,2)/10)) == 0
            toc
            disp(['finishing ',num2str(ceil(iPatt/size(CC.PixelIdxList,2)*100)),'%'])
        end
        if boundPatt(iPatt,6)<thresholdTime
            continue
        end
        if boundPatt(iPatt,4)<thresholdSize | ...
                boundPatt(iPatt,5)<thresholdSize
            continue
        end
        
        durationValid(patCount) = boundPatt(iPatt,6) ;
        sizeValid(patCount) = areaPatt(iPatt) ;
        pattStValid(patCount) = boundPatt(iPatt,3) +0.5 ;
        
        tempBinary = zeros(size(sigIn)) ;
        tempBinary(CC.PixelIdxList{iPatt}) = 1 ;
        
        timeCount = 1 ;
        for iTime = boundPatt(iPatt,3)+0.5:...
                boundPatt(iPatt,3)+boundPatt(iPatt,6)
            sig2D = tempBinary(:,:,iTime) ;
            intCC = bwconncomp(sig2D) ;
            % number of inst. patterns
            intNumPatt{patCount}(timeCount,:) = length(intCC.PixelIdxList) ;
            % inst. size of patterns
            intSizeTotal{patCount}(timeCount,:) = sum(sum(sig2D)) ;
            % inst. centre of patterns
            S = regionprops(sig2D,sigIn(:,:,iTime),{'Centroid','WeightedCentroid'} );
            Centroids{patCount}(timeCount,:) = cat(1, S.Centroid);
            WCentroids{patCount}(timeCount,:) = cat(1, S.WeightedCentroid) ;
            timeCount = timeCount+1 ;
        end 
        patCount = patCount+1 ;
    end
    
    toc
    end
    