% Figure S3B: Changes of pulsing descriptors
addpath('./Functions/')

sites_all = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:69];
nCelltype = 6;
possible_doses = [0 2.5 5 10 20 50 100];

puls_thres = 0.3;

[parentdir,~,~]=fileparts(pwd);
load(fullfile(parentdir,'rawdata','Workspaces','scores_puls_corrected_retracked_all_cleaned_newBTC'));

ylabels = {'','','log_{10} Cyt/Nuc','min'};

resort = [4 1 nan 2 3 6 5]; % Relative to platemap

medians = nan(length(possible_doses),nCelltype,length(features)+1);
highdoses = [];
for isite = sites_all
    sprop = siteprop(isite);
    doseind = sprop.lig_dose == possible_doses;
    if sprop.lig_dose == 100
        highdoses = [highdoses isite];
    end
    
    medians(doseind,resort(sprop.lig_index),1:7) = nanmean(scores_puls(celltypes == isite & scores_puls(:,1)' > puls_thres,:),1);
    medians(doseind,resort(sprop.lig_index),8) = sum(scores_puls(celltypes == isite,1) > puls_thres)/sum(celltypes == isite);
end
medians(1,5,:) = nanmean(medians(1,:,:),2);

resort = [2 3 4 1 6 5];
highdoses = highdoses(resort);

figure

myfeat = [2 3 4 6];

rowstocols = 0.5;
nrows = ceil(length(myfeat)^rowstocols);
ncols = ceil(length(myfeat) / nrows);

colmap = [linspace(0,1,length(highdoses)+1)' ones(length(highdoses)+1,1) ones(length(highdoses)+1,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
markers = {'o','s','v','d','^','>'};
resort2 = [6 2 3 4 5 1];
myslopes_feat = nan(length(highdoses),max(myfeat));
mystd_feat = nan(length(highdoses),max(myfeat));
for iplot = 1:length(myfeat)
    icol = myfeat(iplot);

    subplot(nrows,ncols,iplot)
    hold on
    mycolor = lines(nrows);
    legh = [];
    legstr = cell(length(highdoses),1);
    for irow = 1:length(highdoses)
        isite = highdoses(resort2(irow));
        sprop = siteprop(isite);
        legstr{isite == highdoses} = sprop.lig_name(1:3);
        myind = ~isnan(medians(:,irow,icol));
        legh = [legh plot(0:sum(myind)-1,medians(myind,irow,icol),markers{irow},'MarkerFaceColor',colmap(irow,:),'MarkerEdgeColor',colmap(irow,:),'MarkerSize',6)];

        [axb,s] = polyfit(0:sum(myind)-1,medians(myind,irow,icol)',1);
        plot(0:sum(myind)-1,(0:sum(myind)-1)*axb(1) + axb(2),'-','Color',colmap(irow,:),'LineWidth',2);

        plot([0 sum(myind)-1],[nanmean(scores_puls(scores_puls(:,1)' <= puls_thres,icol),1) nanmean(scores_puls(scores_puls(:,1)' <= puls_thres,icol),1)],'k--','LineWidth',2)
        myslopes_feat(irow,icol) = axb(1);
        Rinv = inv(s.R);
        covmat = sqrt((Rinv*Rinv')*s.normr^2/(s.df));
        mystd_feat(irow,icol) = covmat(1);
    end
    title(features{icol} )
    set(gca,'XLim',[-.5 length(possible_doses)-.5])
    set(gca,'XTick',0:length(possible_doses)-1,'XTickLabel',possible_doses)
    xlabel('Ligand dose [ng/ml]')
    ylabel(ylabels{iplot})
    
end

h = legend(legh,legstr);
ch = get(h,'child');
for ileg = 1:length(ch)/3
    ilegch = (ileg-1)*3+2;
    set(ch(ilegch),'LineStyle','-','LineWidth',2,'Color',colmap(size(colmap,1)-ileg+1,:)); 
end
