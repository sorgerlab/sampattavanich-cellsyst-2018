% Figure 5C: F3aN400 pulsing
addpath('./Functions/')
[parentdir,~,~]=fileparts(pwd);
load(fullfile(parentdir,'rawdata','Workspaces','dists_04182014'));

extension = '04-18-2014';
sites_all = [1:39 41:72];

nrows = 6;
ncols = 12;

puls_thres = .3;

ratio_puls = nan(size(sites_all));
ratio_puls_mat = nan(ncols,nrows);

for isite = sites_all
    ratio_puls(isite) = sum(dists(celltype == isite) > puls_thres) / sum(celltype == isite);
    ratio_puls_mat(subplotpos(isite,ncols)) = ratio_puls(isite);
end

[ratio_puls_sorted,ind_puls_sorted] = sort(ratio_puls,'descend');

[X,Y] = meshgrid(1:size(ratio_puls_mat,2), 1:1:size(ratio_puls_mat,1)); 

valid = ~isnan(ratio_puls_mat); 
M = griddata(X(valid),Y(valid),ratio_puls_mat(valid),X,Y);
zlab = 'Ratio of pulsing cells';

colmap = cbrewer('seq','YlGnBu',201);

figure

imagesc(M(1:6,:))
colormap(colmap)
set(gca,'XDir','Reverse')
title(['MCF10A - ' zlab])
xlabel('EGF [ng/mL]')
ylabel('PD0325901 [\muM]')
set(gca,'XTick',1:6,'XTickLabel',[100 20 4 .8 .16 0])
set(gca,'YTick',1:6,'YTickLabel',[.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0]*1000)
colorbar

figure

imagesc(M(12:-1:7,:))
colormap(colmap)
set(gca,'XDir','Reverse')
title(['184A1 - ' zlab])
xlabel('EGF [ng/mL]')
ylabel('PD0325901 [\muM]')
set(gca,'XTick',1:6,'XTickLabel',[100 20 4 .8 .16 0])
set(gca,'YTick',1:6,'YTickLabel',[.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0]*1000)
colorbar
