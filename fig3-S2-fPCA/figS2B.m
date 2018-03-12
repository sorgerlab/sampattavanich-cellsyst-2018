% Figure S2B: Pair-wise comparison of fPCA scores

[parentdir,~,~]=fileparts(pwd);
load(fullfile(parentdir,'rawdata','Workspaces','scores_early_5basis_noFGF_newBTC'));

sites_all = [17:-1:11; 37:-1:31; 44:50; 4:10; 64:69 10; 57:-1:51];
sites_high = [4 17 37 44 57 64];
sites_unst = [10 11 31 50 51];

myind = ~isnan(scores_early(1,:));

ps = 1:3;

colmap = hsv(size(sites_all,1));

for ipc = ps
    pvals_all = nan(size(sites_all));
    for i = 1:length(sites_all(:))
        isite = sites_all(i);
        pvals_all(i) = ranksum(scores_early(ipc,ismember(celltypes,sites_unst) & myind),scores_early(ipc,celltypes == isite & myind));
    end

    figure
    hold on
    legh = [];
    legstr = {};
    for isite = 1:size(sites_all,1)
        legh = [legh plot(pvals_all(isite,:),'o-','Color',colmap(isite,:))];
        s = siteprop(sites_all(isite,1));
        legstr{end+1} = s.lig_name;
    end
    legend(legh,legstr,'Location','SouthWest')
    plot(get(gca,'XLim'),[1e-10 1e-10],'k--')
    if ipc == 1
        set(gca,'YLim',[1e-15 1])
    end
    set(gca,'YScale','log','XDir','reverse')
    title(sprintf('PC %i',ipc))
    xlabel('ligand dose')
    ylabel('p-value')
    set(gca,'XTick',1:7,'XTickLabel',[100 50 20 10 5 2.5 0])
end
