function [gm,labeling,surf_lh] = genGradients
% generate gradients
addpath(genpath([pwd,'/ToolOthers/fMRI/BrainSpace-master/matlab']));

% First load mean connectivity matrix and Schaefer parcellation
conn_matrix = load_group_fc('schaefer',400);
labeling = load_parcellation('schaefer',400);

% The loader functions output data in a struct array for when you
% load multiple parcellations. Let's just bring them to numeric arrays.
conn_matrix = conn_matrix.schaefer_400;
labeling = labeling.schaefer_400;

% and load the conte69 hemisphere surfaces
[surf_lh, surf_rh] = load_conte69();


% h = plot_hemispheres(labeling, {surf_lh,surf_rh});
% colormap(h.figure,lines(401))



% Construct the gradients
gm = GradientMaps();
gm = gm.fit(conn_matrix);



% figure;plot_hemispheres(gm.gradients{1}(:,1:2),{surf_lh,surf_rh}, ...
%              'parcellation', labeling, ...
%              'labeltext',{'Gradient 1','Gradient 2'});
%          
%          figure;plot_hemispheres(gm.gradients{1}(1:size(gm.gradients{1},1)/2,1),{surf_lh}, ...
%              'parcellation', labeling(1:size(labeling,1)/2), ...
%              'labeltext',{'Gradient 1'});
