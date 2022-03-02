%% visualize amplitude patterns of different bands

fileDir = [pwd,'/preprocess/Preprocess_BP/'] ;

fileName1 = dir(fileDir) ;
iSub = 1 ;
fileNum = 2 + iSub ;  % find the resting data

load([fileDir,fileName1(fileNum).name])
sigBPass_bp = sigBPass ;
params_bp = params ;

%%
fileDir = [pwd,'/preprocess/Preprocess_eigenWide/'] ;
% fileDir = [pwd,'/preprocess/Preprocess_eigen/'] ;


fileName1 = dir(fileDir) ;
iSub = 1 ;
fileNum = 2 + iSub ;  % find the resting data

load([fileDir,fileName1(fileNum).name])
sigBPass_eg = sigBPass ;
params_eg = params ;

%% visualize

flagPlot = 1 ;
if flagPlot
flagSaveVideo = 1 ;
if flagSaveVideo
    t =  datetime('now') ;
    dateStr = datestr(t,'mmmmdd_HH:MM') ;
    vidTitle = [pwd,'/Results/Project3/movies/'] ;
    saveFileName = [vidTitle,'fMRI_s1','amp_compare_eigenWide_',dateStr] ;
    vidObj = VideoWriter(saveFileName,'Motion JPEG AVI');
    v.Quality = 50 ;
    vidObj.FrameRate = 5 ;  % 5
    open(vidObj);
    h.figure = figure('Color','white') ;
    % set(gcf,'Position', [224 80 848 604])
    % set(gcf,'Position',[358 266 560 420])
    set(gcf,'Position',[1 50 1280 636]) % 
else
    h.figure = figure('Color','white');
    set(gcf,'Position', [224 80 848 604])
    set(gcf,'Position',[1 50 1280 636])
end
cScale = 0.4 ;
cmin1 = cScale*min(sigBPass_bp(1,:,:,:),[],'all') ;
cmax1 = cScale*max(sigBPass_bp(1,:,:,:),[],'all') ;
cmax1 = -cScale*min(sigBPass_bp(1,:,:,:),[],'all') ;

cmin2 = cScale*min(sigBPass_bp(2,:,:,:),[],'all') ;
cmax2 = cScale*max(sigBPass_bp(2,:,:,:),[],'all') ;
cmax2 = -cScale*min(sigBPass_bp(2,:,:,:),[],'all') ;

cmin3 = cScale*min(sigBPass_bp(3,:,:,:),[],'all') ;
cmax3 = cScale*max(sigBPass_bp(3,:,:,:),[],'all') ;
cmax3 = -cScale*min(sigBPass_bp(3,:,:,:),[],'all') ;

cmin4 = cScale*min(sigBPass_eg(1,:,:,:),[],'all') ;
cmax4 = cScale*max(sigBPass_eg(1,:,:,:),[],'all') ;
cmax4 = -cScale*min(sigBPass_eg(1,:,:,:),[],'all') ;

cmin5 = cScale*min(sigBPass_eg(2,:,:,:),[],'all') ;
cmax5 = cScale*max(sigBPass_eg(2,:,:,:),[],'all') ;
cmax5 = -cScale*min(sigBPass_eg(2,:,:,:),[],'all') ;

cmin6 = cScale*min(sigBPass_eg(3,:,:,:),[],'all') ;
cmax6 = cScale*max(sigBPass_eg(3,:,:,:),[],'all') ;
cmax6 = -cScale*min(sigBPass_eg(3,:,:,:),[],'all') ;


vDownR = 1 ;
cb = [0.8,0.8,0.8] ;
for iTime = 1:1200
    subplot(2,3,1)
    imagesc(squeeze(sigBPass_bp(1,:,:,iTime)))
        caxis([cmin1,cmax1])
    set(gca,'yDir','normal')
    hold on
    plot22SecBound(cb)
    title(['Band 1, sigma 64-32 time = ',num2str(iTime)])
    colorbar

    
    subplot(2,3,2)
    imagesc(squeeze(sigBPass_bp(2,:,:,iTime)))
    caxis([cmin2,cmax2])
    set(gca,'yDir','normal')
    hold on
    plot22SecBound(cb)
    title(['Band 2, sigma 32-16'])
    colorbar

    subplot(2,3,3)
    imagesc(squeeze(sigBPass_bp(3,:,:,iTime)))
    caxis([cmin3,cmax3])
    set(gca,'yDir','normal')
    hold on
    plot22SecBound(cb)
    title(['Band 3, sigma 16-8'])
    colorbar
    
    subplot(2,3,4)
    imagesc(squeeze(sigBPass_eg(1,:,:,iTime)))
    caxis([cmin4,cmax4])
    set(gca,'yDir','normal')
    hold on
    plot22SecBound(cb)
    title(['Eigen 1&2, sigma 84 - 29'])
    colorbar
    
    subplot(2,3,5)
    imagesc(squeeze(sigBPass_eg(2,:,:,iTime)))
    caxis([cmin5,cmax5])
    set(gca,'yDir','normal')
    hold on
    plot22SecBound(cb)
    title(['Eigen 3&4, sigma 29 - 18'])
    colorbar
    
    subplot(2,3,6)
    imagesc(squeeze(sigBPass_eg(3,:,:,iTime)))
    caxis([cmin6,cmax6])
    set(gca,'yDir','normal')
    hold on
    plot22SecBound(cb)
    title(['Eigen 5&6, sigma 18 - 13'])
    colorbar
    
    if flagSaveVideo
        writeVideo(vidObj, im2frame(print(gcf,'-RGBImage')));
    else
        pause(0.1)
    end
    clf
end
if flagSaveVideo
    close(vidObj)
end
end