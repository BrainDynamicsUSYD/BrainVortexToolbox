function sigSSpec_reshape=spaceFreq_fMRI(posValid,dataIn,xCord,yCord,...
    flagPlotSpe, flagCheckIntp,mask)

%% re-interpolate data (time consuming)
[X,Y] = meshgrid(xCord, yCord);
sigSSpec = dataIn ;     

sigSSpec_reshape = [] ;
disp('start interpolation (downsample)...')
tic
if sum(mask(:))~= 0
    flagMask = 1 ;
end
for iTime=1:size(sigSSpec,2)    
%      sigSSpec_reshape(:,:,iTime) = griddata(double(posValid.vertices(1:10:end,1)),...
%         double(posValid.vertices(1:10:end,2)),sigSSpec(1:10:end,iTime),X,Y,'cubic') ;
        sigSSpec_reshape(:,:,iTime) = griddata(double(posValid.vertices(:,1)),...
            double(posValid.vertices(:,2)),sigSSpec(:,iTime),X,Y,'cubic') ;
    if flagMask
        sigSSpec_reshape(:,:,iTime) = sigSSpec_reshape(:,:,iTime).*mask ;
    end
end
toc

%% plot the spatial spectrum 
if flagPlotSpe
[Pxx,F] = pmtm(zscore(sigSSpec_reshape(100,:,1)),[],[],1e3) ;
figure;
coeffs = polyfit(log(F(F>10)),log(Pxx(F>10)),1); 
fittedX = log(F);
fittedLogY = polyval(coeffs, fittedX);
fittedY = exp(fittedLogY) ;

loglog(F,Pxx)
hold on
loglog(F,fittedY)
end

%% check the re-interpolation by visualisation
if flagCheckIntp
    %% original data
sigIn = sigSSpec ;  % angle(sigSSpec) ;
numColor = 1000 ;
tempColor = jet(numColor) ;
maxData = max(sigIn(:)) ;
minData = min(sigIn(:)) ;
colorAxis = linspace(minData,maxData/2,numColor) ;
idxColor1 = zeros(size(sigIn)) ;
disp('start visulisation 1 ...')
tic
for iTime = 1:size(sigIn,2)
    % [cordI,cordJ] = ind2sub(size(sigIn),iData) ;
    [~,idxColor1(:,iTime)] = min(abs(sigIn(:,iTime)-ones(size(sigIn,1),1)*colorAxis),[],2) ;
    % disp([num2str(iTime/size(sigIn,2)),'%'])
end
toc
    %% interpolated data
sigIn = sigSSpec_reshape ;  % angle(sigSSpec_reshape) ;  % notice the size change
% numColor = 1000 ;
% tempColor = jet(numColor) ;
% maxData = max(sigIn(:)) ;
% minData = min(sigIn(:)) ;
% colorAxis = linspace(minData,maxData/2,numColor) ;
sigIn = reshape(sigIn,size(sigIn,1)*size(sigIn,2),size(sigIn,3)) ;
idxColor2 = zeros(size(sigIn)) ;
disp('start visulisation 2 ...')
tic
for iTime = 1:size(sigIn,2)
    % [cordI,cordJ] = ind2sub(size(sigIn),iData) ;
    [~,idxColor2(:,iTime)] = min(abs(sigIn(:,iTime)-ones(size(sigIn,1),1)*colorAxis),[],2) ;
    % disp([num2str(iTime/size(sigIn,2)),'%'])
end
toc
X_list = reshape(X,size(X,1)*size(X,2),1) ;
Y_list = reshape(Y,size(Y,1)*size(Y,2),1) ;

iTime = 100 ;
figure;
scatter(posValid.vertices(:,1),posValid.vertices(:,2),4,tempColor(idxColor1(:,iTime),:)) ;
xlim([xCord(1) xCord(end)])
ylim([yCord(1) yCord(end)])

figure;
scatter(X_list,Y_list,4,tempColor(idxColor2(:,iTime),:)) ;


end
