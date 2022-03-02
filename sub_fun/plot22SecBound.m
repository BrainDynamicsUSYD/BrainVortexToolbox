function plot22SecBound(colorPlot)

addpath(genpath([pwd,'/Toolbox_CSC']))
addpath(genpath([pwd,'/ToolOthers/nanconv']))
% set basic parameters
params.sigmScale = 2.^(6:-1:2) ;              % log spacing of scales
params.downSRate = 2 ;                        % downsample the re-interpolation
params.xCord = -250:params.downSRate:250 ;    % coordinate re-interpolation
params.yCord = -150:params.downSRate:200 ;
params.fsTem = 1/0.72 ;                       % temporal sampling rate
params.zscore = 1 ;                           % flag, 1 for zscore, 0 for no     


% load HCP
addpath(genpath([pwd,'/ToolOthers/fMRI/cifti-matlab-master']))

folderName = [pwd,'/Data/fMRI2/100610_3T_structural/MNINonLinear/fsaverage_LR32k/'] ;
fileName = '100610.L.flat.32k_fs_LR.surf.gii' ;
posFile = [folderName,fileName] ;
posLC = gifti(posFile) ;

parLabel2 = ft_read_cifti([pwd,'/Data/HCP/ResultsRegions_ROI.dlabel.nii']) ;
parValid = parLabel2.x1(1:32492) ;
parValid(isnan(parLabel2.x1(1:32492)),:) = [] ;
posValid = posLC ;
posValid.vertices(isnan(parLabel2.x1(1:32492)) ,:) = [] ; 

hold on;
for iPar = setdiff(1:22, [12,13,18,20])
    idxA4 = find(parValid==iPar);
    
    y = (double(posValid.vertices(idxA4,1)) - min(params.xCord)+1)/params.downSRate+0.5 ;
    x = (double(posValid.vertices(idxA4,2)) - min(params.yCord)+1)/params.downSRate+0.5 ;
    k = boundary(x,y) ;
    hold on;
    lh = plot(  (y(k))+0.5,(x(k))+0.5,'--','color',colorPlot,'lineWidth',2);
    % lh.Color = [lh.Color 0.4] ;
end 

end