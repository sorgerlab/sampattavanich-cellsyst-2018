% Figure 3A: Decomposing early trajectories by fPCA
addpath('./Functions/')
close all
[parentdir,~,~]=fileparts(pwd);
load(fullfile(parentdir,'rawdata','Workspaces',['site_' num2str(isite)]));
load(fullfile(parentdir,'rawdata','Workspaces','harm_basis_fPCA_5basis_noFGF_newBTC_rot'));
load(fullfile(parentdir,'rawdata','Workspaces','scores_early_5basis_noFGF_newBTC'));
myextension = '130722_corrected_retracked_all_cleaned';

propvar = squeeze(sum(harmscr.^2));
propvar = propvar./sum(propvar);
propvar = propvar.*sum(c_signal_pcastr.varprop(1:length(propvar)));
propvar = round(propvar*1e4)./1e2;

sites_all = [17 57 64];
sigs = [26 35 99];
colind = [1 6 5];
basisind = 1:3; % Only use basis 1-3 to fit data

legstr = cell(1,length(sites_all));
c_signal_single = [];
scores_single = nan(size(scores_early,1),length(sigs));

for icount = 1:length(sites_all)
    isite = sites_all(icount);
    load(fullfile(parentdir,'rawdata','Workspaces',['site_' num2str(isite) '_' myextension]));
    c_signal_single(:,icount) = log10(intensity(:,sigs(icount)));
    scores_tmp = scores_early(:,celltypes==isite);
    scores_single(:,icount) = scores_tmp(:,sigs(icount));
    s = siteprop(isite);
    legstr{icount} = s.lig_name;
end

nrows = 3;
ncols = 2;

colmap = [linspace(0,1,7)' ones(7,1) ones(7,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
markers = {'o','s','v','d','^','>'};

time_range = [50.7 197.8];
times_fine = linspace(time_range(1),time_range(2),501);

basis_eval = eval_basis(harm_basis,times_fine);
times_fine(end+1) = 200;
basis_eval(end+1,:) = basis_eval(end,:);

figure
hold on
legh = [];
for iplot = basisind
    legh = [legh plot(timestamp,c_signal_single(:,iplot),markers{colind(iplot)},'Color',colmap(colind(iplot),:))];
    plot(timestamp,timestamp*0,'k:')
    plot(times_fine,basis_eval(:,basisind)*scores_single(basisind,iplot),'Color',colmap(colind(iplot),:))
end
set(gca,'XLim',[50 200],'XTick',60:30:180,'XTickLabel',-60:30:60)

h = legend(legh,legstr);
ch = get(h,'child');
for ileg = 1:length(ch)/3
    ilegch = (ileg-1)*3+2;
    set(ch(ilegch),'LineStyle','-','Color',colmap(colind(end-ileg+1),:)); 
end

figure

for iplot = basisind
    subplot(nrows,ncols,ncols*(iplot-1)+ncols-1)
    ylabel(['Harmonic ' num2str(iplot)])
    xlabel('time [min]')
    hold on

    plot(times_fine,basis_eval(:,iplot),'k')   
    plot(time_range,[0 0],'k:')
    
    set(gca,'XLim',[50 200],'XTick',60:30:180,'XTickLabel',-60:30:60)
    set(gca,'YLim',[-.22 .22])
    
    text(60,-.14,sprintf('Explains %g%% variance',propvar(iplot)))
end

ylims = [[-.2 .1];[-.1 .15];[-.05 .15]];
for iplot = basisind
    subplot(nrows,ncols,ncols*(iplot-1)+ncols)
    ylabel('Score')
    hold on
    
    for ilig = 1:size(scores_single,2)
        h = bar(ilig,scores_single(iplot,ilig));
        set(h,'FaceColor',colmap(colind(ilig),:))
    end
    
    set(gca,'XTick',1:size(scores_single,2),'XTickLabel',legstr)
    set(gca,'YLim',ylims(iplot,:))
end
