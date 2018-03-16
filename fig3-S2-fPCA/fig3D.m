% Figure 3D: Medians of fPC scores
addpath('./Functions/')

sites_all = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:70];
nCelltype = 6;
possible_doses = [0 2.5 5 10 20 50 100];

[parentdir,~,~]=fileparts(pwd);
load(fullfile(parentdir,'rawdata','Workspaces','scores_early_5basis_noFGF_newBTC'));

resort = [4 1 nan 2 3 6 5]; % Relative to platemap

medians = nan(length(possible_doses),nCelltype,5);
highdoses = [];
for isite = sites_all
    sprop = siteprop(isite);
    doseind = sprop.lig_dose == possible_doses;
    if sprop.lig_dose == 100
        highdoses = [highdoses isite];
    end
    
    medians(doseind,resort(sprop.lig_index),:) = nanmedian(scores_early(:,celltypes == isite),2);
end
medians(1,5,:) = nanmean(medians(1,:,:),2);

resort = [2 3 4 1 6 5];
highdoses = highdoses(resort);

nrows = 1;
ncols = 3;

figure

setFigure(gcf,1,.5,12)

subplot(nrows,ncols,1)

hold on
icol = 1;
legh = [];
colmap = [linspace(0,1,length(highdoses)+1)' ones(length(highdoses)+1,1) ones(length(highdoses)+1,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
markers = {'o','s','v','d','^','>'};
for irow = 1:length(highdoses)
    plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6)

    [axb,s] = polyfit(0:length(possible_doses)-1,medians(1:end,irow,icol)',1);
    plot(0:length(possible_doses)-1,(0:length(possible_doses)-1)*axb(1) + axb(2),'-','Color',colmap(irow,:),'LineWidth',2);
end
title('Harmonic 1')
set(gca,'XLim',[-.5 length(possible_doses)-.5])
set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
xlabel('Ligand dose [ng/ml]')


subplot(nrows,ncols,2)
hold on
icol = 2;
legh = [];
for irow = 1:length(highdoses)
    plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6)

    [axb,s] = polyfitZero(1:length(possible_doses)-1,medians(2:end,irow,icol)'-mean(medians(1,:,icol)),1);
    legh = [legh plot(0:length(possible_doses)-1,(0:length(possible_doses)-1)*axb(1) + mean(medians(1,:,icol)),'-','Color',colmap(irow,:),'LineWidth',2)];
end
title('Harmonic 2')
set(gca,'XLim',[-.5 length(possible_doses)-.5])
set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
xlabel('Ligand dose [ng/ml]')


s5 = subplot(nrows,ncols,3);
hold on
icol = 3;
mycolor = lines(nrows);
legh = [];
legstr = cell(length(highdoses),1);
resort2 = [6 2 3 4 5 1];
for irow = 1:length(highdoses)
    isite = highdoses(resort2(irow));
    sprop = siteprop(isite);
    legstr{isite == highdoses} = sprop.lig_name(1:3);
    legh(irow) = plot(0:length(possible_doses)-1,medians(:,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6);

    [axb,s] = polyfitZero(1:length(possible_doses)-1,medians(2:end,irow,icol)'-mean(medians(1,:,icol)),1);
    plot(0:length(possible_doses)-1,(0:length(possible_doses)-1)*axb(1) + mean(medians(1,:,icol)),'-','Color',colmap(irow,:),'LineWidth',2);
end

h = legend(legh,legstr,'Location','NorthWest');
ch = get(h,'child');
for ileg = 1:length(ch)/3
    ilegch = (ileg-1)*3+2;
    set(ch(ilegch),'LineStyle','-','LineWidth',2,'Color',colmap(size(colmap,1)-ileg+1,:)); 
end

title('Harmonic 3')
set(gca,'XLim',[-.5 length(possible_doses)-.5])
set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
xlabel('Ligand dose [ng/ml]')
