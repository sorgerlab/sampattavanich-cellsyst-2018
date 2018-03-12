% Figures 5B and S4C: Pharmacological inhibition
addpath('./Functions/')
[parentdir,~,~]=fileparts(pwd);
load(fullfile(parentdir,'rawdata','Workspaces','scores_early_5basis_noFGF_AKTi'));
meki = load(fullfile(parentdir,'rawdata','Workspaces','scores_early_5basis_noFGF_MEKi'));
mekipuls = load(fullfile(parentdir,'rawdata','Workspaces','scores_puls_corrected_retracked_all_cleaned_newBTC_MEKi'));
load(fullfile(parentdir,'rawdata','Workspaces','scores_puls_corrected_retracked_all_cleaned_newBTC_ATKi'));
noInh = load(fullfile(parentdir,'rawdata','Workspaces','scores_early_5basis_noFGF_newBTC'));
noInhpuls = load(fullfile(parentdir,'rawdata','Workspaces','scores_puls_corrected_retracked_all_cleaned_newBTC'));

pcs = [2 3];

sites_all = [37 44 4 64];
sites_akti = [39 42 2 62];
sites_meki = [40 41 1 61];

puls_thres = .3;

highinds = ismember(celltypes,sites_akti);
scores_puls = scores_puls(highinds,1);
highindsNoInh = ismember(noInhpuls.celltypes,sites_all);
noInhpuls.scores_puls = noInhpuls.scores_puls(highindsNoInh,1);
highindsMEKi = ismember(mekipuls.celltypes,sites_meki);
mekipuls.scores_puls = mekipuls.scores_puls(highindsMEKi,1);

figure

plot(noInh.scores_early(pcs(1),highindsNoInh),noInh.scores_early(pcs(2),highindsNoInh),'.','Color',[.7 .7 .7]) % WT
hold on
plot(meki.scores_early(pcs(1),highindsMEKi),meki.scores_early(pcs(2),highindsMEKi),'.','Color',[.7 .7 .7]) % MEKi

colmap = hsv(5+1);
markers = {'o','s','v','d','^','>'};
colmap = colmap(2:end,:);
markers = markers(2:end);
legstr = {};
for i = 1:length(sites_akti)
    isite = sites_akti(i);
    s = siteprop(isite);
    titstr = s.lig_name;
    legstr{end+1} = titstr;
    
    plot(nanmean(scores_early(pcs(1),celltypes == isite)),nanmean(scores_early(pcs(2),celltypes == isite)),markers{isite == sites_akti},'Color',colmap(isite == sites_akti,:),'MarkerFaceColor',colmap(isite == sites_akti,:),'MarkerEdgeColor','w','MarkerSize',12)
    
    isite3 = sites_meki(i);
    plot(nanmean(meki.scores_early(pcs(1),meki.celltypes == isite3)),nanmean(meki.scores_early(pcs(2),meki.celltypes == isite3)),markers{isite3 == sites_meki},'Color',colmap(isite3 == sites_meki,:),'MarkerFaceColor',colmap(isite3 == sites_meki,:),'MarkerEdgeColor','k','MarkerSize',12)
    isite2 = sites_all(i);
    plot(nanmean(noInh.scores_early(pcs(1),noInh.celltypes == isite2)),nanmean(noInh.scores_early(pcs(2),noInh.celltypes == isite2)),markers{isite2 == sites_all},'Color',colmap(isite2 == sites_all,:),'MarkerFaceColor',colmap(isite2 == sites_all,:),'MarkerSize',12)
    
    plot([nanmean(noInh.scores_early(pcs(1),noInh.celltypes == isite2)) nanmean(meki.scores_early(pcs(1),meki.celltypes == isite3))],[nanmean(noInh.scores_early(pcs(2),noInh.celltypes == isite2)) nanmean(meki.scores_early(pcs(2),meki.celltypes == isite3))],'-','Color',colmap(isite == sites_akti,:),'LineWidth',2)
    plot([nanmean(scores_early(pcs(1),celltypes == isite)) nanmean(noInh.scores_early(pcs(1),noInh.celltypes == isite2))],[nanmean(scores_early(pcs(2),celltypes == isite)) nanmean(noInh.scores_early(pcs(2),noInh.celltypes == isite2))],'k--','Color','k')
end
set(gca,'XLim',[-0.3 0.25],'YLim',[-0.02 0.1])

ylabel(['PC ' num2str(pcs(2))])
xlabel(['PC ' num2str(pcs(1))])

set(gca,'CLim',[0 1])
colormap(colmap(1:end-1,:))
colorbar('YTick',linspace(1./(2*length(sites_akti)),1-1./(2*length(sites_akti)),length(sites_akti)),'YTickLabel',legstr,'TickLength', 0) % Vertical colorbar

figure
hold on

rng(0) % Make sure that bootstrap samples are reproducible
ratio_fun = @(x) sum(x > puls_thres) / length(x);
errorb = nan(3*length(sites_akti),4);

legh = [];
legstr = {};
for i = 1:length(sites_akti)
    isite = sites_akti(i);
    s = siteprop(isite);
    
    iinh = 3;
    errorb((i-1)*3+iinh,1) = mod(iinh-1,3)*(length(sites_akti)+1)+i;
    errorb((i-1)*3+iinh,2) = sum(scores_puls(celltypes(highinds) == isite) > puls_thres)./sum(celltypes(highinds) == isite);
    errorb((i-1)*3+iinh,3:4) = bootci(2000,{ratio_fun,scores_puls(celltypes(highinds) == isite)},'alpha',.32);
    legh = [legh bar(errorb((i-1)*3+iinh,1),errorb((i-1)*3+iinh,2),'FaceColor',colmap(isite == sites_akti,:))];
    iinh = 1;
    isiteH = sites_all(i);
    errorb((i-1)*3+iinh,1) = mod(iinh-1,3)*(length(sites_akti)+1)+i;
    errorb((i-1)*3+iinh,2) = sum(noInhpuls.scores_puls(noInhpuls.celltypes(highindsNoInh) == isiteH) > puls_thres)./sum(noInhpuls.celltypes(highindsNoInh) == isiteH);
    errorb((i-1)*3+iinh,3:4) = bootci(2000,{ratio_fun,noInhpuls.scores_puls(noInhpuls.celltypes(highindsNoInh) == isiteH)},'alpha',.32);
    bar(errorb((i-1)*3+iinh,1),errorb((i-1)*3+iinh,2),'FaceColor',colmap(isiteH == sites_all,:));
    iinh = 2;
    isiteM = sites_meki(i);
    errorb((i-1)*3+iinh,1) = mod(iinh-1,3)*(length(sites_akti)+1)+i;
    errorb((i-1)*3+iinh,2) = sum(mekipuls.scores_puls(mekipuls.celltypes(highindsMEKi) == isiteM) > puls_thres)./sum(mekipuls.celltypes(highindsMEKi) == isiteM);
    errorb((i-1)*3+iinh,3:4) = bootci(2000,{ratio_fun,mekipuls.scores_puls(mekipuls.celltypes(highindsMEKi) == isiteM)},'alpha',.32);
    bar(errorb((i-1)*3+iinh,1),errorb((i-1)*3+iinh,2),'FaceColor',colmap(isiteM == sites_meki,:));
    legstr{end+1} = s.lig_name;
end

title('184A1')
set(gca,'XTick',3:5:13,'XTickLabel',{'WT','MEKi','AKTi'});

legend(legh,legstr)

ylabel('fraction of pulsing cells')

errorbar(errorb(:,1),errorb(:,2),errorb(:,2)-errorb(:,3),errorb(:,4)-errorb(:,2),'LineStyle','none','Color','k');

sites_single = [sites_all sites_akti sites_meki];
celltypes_single = [];
c_signal_single = [];

for isite = sites_single
    tmp = load(fullfile(parentdir,'rawdata','Workspaces',sprintf('site_%d_130722_corrected_retracked_all_paper_cleaned',isite)));
    c_signal_single = [c_signal_single tmp.intensity];
    timestamp_single = tmp.timestamp;
    celltypes_single = [celltypes_single ones(1,size(tmp.intensity,2))*isite];
end

nrows = 2;
ncols = 3;

ntraces = 25; % This is the maximal number of traces to be plotted

figure
colmap = hsv(5+1);
colmap = colmap(2:end,:);
for i = 1:length(sites_all)
    subplot(nrows,ncols,i)
    hold on

    ex2 = find(celltypes_single==sites_akti(i) & sum(isnan(c_signal_single),1) < 10);
    plot(timestamp_single-120,c_signal_single(:,ex2(1:min([ntraces length(ex2)]))),'Color',[.7 .7 .7])
    
    ex1 = find(celltypes_single == sites_all(i) & sum(isnan(c_signal_single),1) < 10);
    plot(timestamp_single-120,c_signal_single(:,ex1(1:min([ntraces length(ex1)]))),'Color',colmap(i,:))

    a = plot(timestamp_single-120,nanmean(c_signal_single(:,celltypes_single==sites_all(i)),2),'k','LineWidth',2);
    b = plot(timestamp_single-120,nanmean(c_signal_single(:,celltypes_single==sites_akti(i)),2),'k--','LineWidth',2);
    
    s = siteprop(sites_all(i));
    titstr = s.lig_name;
    title(titstr)
    
    legend([a b],{'Ligand','+ AKTi'},'Location','SouthEast')
    
    set(gca,'YLim',[.9 1.1])
    set(gca,'XLim',[-70 720])
end

figure
for i = 1:length(sites_all)
    subplot(nrows,ncols,i)
    hold on

    ex2 = find(celltypes_single==sites_meki(i) & sum(isnan(c_signal_single),1) < 10);
    plot(timestamp_single-120,c_signal_single(:,ex2(1:min([ntraces length(ex2)]))),'Color',[.7 .7 .7])
    
    ex1 = find(celltypes_single == sites_all(i) & sum(isnan(c_signal_single),1) < 10);
    plot(timestamp_single-120,c_signal_single(:,ex1(1:min([ntraces length(ex1)]))),'Color',colmap(i,:))

    a = plot(timestamp_single-120,nanmean(c_signal_single(:,celltypes_single==sites_all(i)),2),'k','LineWidth',2);
    b = plot(timestamp_single-120,nanmean(c_signal_single(:,celltypes_single==sites_meki(i)),2),'k--','LineWidth',2);
    
    s = siteprop(sites_all(i));
    titstr = s.lig_name;
    title(titstr)
    
    legend([a b],{'Ligand','+ MEKi'},'Location','SouthEast')
    
    set(gca,'YLim',[.9 1.1])
    set(gca,'XLim',[-70 720])
end