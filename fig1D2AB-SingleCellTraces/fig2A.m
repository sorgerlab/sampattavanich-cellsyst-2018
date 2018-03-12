% Figure 2A: F3aN400-Venus trajectories by ligand
addpath('./Functions/')

sites_all = [17 37 44 64 57 4];

times = cell(0);
signals = cell(0);
celltype = [];
legstr = {};

for isite = sites_all
    [parentdir,~,~]=fileparts(pwd);
    load(fullfile(parentdir,'rawdata','Workspaces',['site_' num2str(isite)]));
    times{end+1} = timestamp;  
    signals{end+1} = log10(intensity);
    celltype = [celltype ones(1,size(intensity,2))*isite];
    s = siteprop(isite);
    legstr{end+1} = s.lig_name;
end

timestamp = times{1};
c_signal = cell2mat(signals);

time_range = [50 605];

[~,range_ind_min] = min(abs(timestamp - time_range(1)));
[~,range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max+1;

figure

ncols = 3;
nrows = 2;

colmap = [linspace(0,1,length(sites_all)+1)' ones(length(sites_all)+1,1) ones(length(sites_all)+1,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));

for iplot = 1:length(sites_all)-1
    subplot(nrows,ncols,iplot+1)
    hold on
    
    isite = sites_all(iplot);
    s = siteprop(isite);

    first_n = 7; % Plot first_n traces colored (2 times)

    tmp_puls_strengths = edge_snr_score_pw_distdur(isite,[],0,1/120,'harm_basis_130722_corrected_retracked_all_cleaned_late',1);
    [~,ind_tmp_pul_str] = sort(tmp_puls_strengths);
    ind_isite = [ind_tmp_pul_str(end:-1:end-first_n+1) ind_tmp_pul_str(round(linspace(1,(length(ind_tmp_pul_str)-first_n),first_n)))];

    c_signal_single = c_signal(:,celltype == isite);
    plot(repmat(timestamp(range_ind),1,size(c_signal_single,2)),c_signal_single(range_ind,:),'g','color',[0.7 0.7 0.7])
    plot(repmat(timestamp(range_ind),1,2*first_n),c_signal_single(range_ind,ind_isite))
    title(s.lig_name)

    if iplot == 5
        xlabel('time [min]')
    end
    if iplot == 1
        ylabel('log_{10} FOXO3a [Cyt/Nuc]')
    end

    ylim = [-1 1]*.05;

    set(gca,'XLim',[50 600],'YLim',ylim);
    set(gca,'XTick',120:100:520,'XTickLabel',0:100:400,'YTick',[-0.04 -0.02 0 0.02 0.04])
end

subplot(nrows,ncols,1)
hold on
legh = [];
for iplot = 1:length(sites_all)
    isite = sites_all(iplot);
    legh = [legh plot(timestamp(range_ind),nanmean(c_signal(range_ind,celltype == isite),2),'Color',colmap(iplot,:))];
end
title('Averaged')
set(gca,'XLim',time_range,'YLim',ylim/2)
set(gca,'XTick',120:100:520,'XTickLabel',0:100:400,'YTick',[-0.04 -0.02 0 0.02 0.04])
legend(legh,legstr)