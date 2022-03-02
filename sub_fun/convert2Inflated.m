function convert2Inflated(Vq,posLC_inflat,colMap,colAxis,viewAngle)



S = convert_surface(posLC_inflat);

% h.figure = figure('Color','white','Units','normalized','Position',[0.1828 0.3630 0.3906 0.3787]);
% h.axes = axes('Position',[-.1+1*.133 .9-.2*1 .2 .2]);
h.trisurf = trisurf(S.tri, ...
            S.coord(1,:), ...
            S.coord(2,:), ...
            S.coord(3,:), ...
            Vq, ...
            'EdgeColor', 'None');
 set(gca                    , ...
    'Visible'           , 'off'         , ...
    'DataAspectRatio'   , [1 1 1]       , ...
    'PlotBoxAspectRatio', [1 1 1]           );
        material dull; lighting phong;

set(gca,'View',viewAngle);
%     h.camlight = camlight();
colormap(colMap)
caxis(colAxis)
    
    
    
