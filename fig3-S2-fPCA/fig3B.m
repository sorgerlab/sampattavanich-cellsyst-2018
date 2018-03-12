% Figure 3B: Scores plot fPC2 vs fPC3
addpath('./functions/')

[parentdir,~,~]=fileparts(pwd);
datapath = fullfile(parentdir,'rawdata','Workspaces');
data = loadcsv('130722_SCfeat',datapath);


sites_all = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:69];

celltypes = data(ismember(data(:,2),sites_all),2)';
scores_early = data(ismember(data(:,2),sites_all),19:21)';

outliers = false(size(celltypes));
outlier_thres = .05;
nNeighbours = 8;
for iscore = 1:length(scores_early)
    sorted_dists = sort(sqrt(sum((repmat(scores_early(2:3,iscore),1,size(scores_early,2))-scores_early(2:3,:)).^2,1)));
    outliers(iscore) = sorted_dists(nNeighbours+1) > outlier_thres;
end


highdoses = [];
for isite = sites_all
    sprop = siteprop(isite);

    if sprop.lig_dose == 100
        highdoses = [highdoses isite];
    end
end

resort = [2 3 4 1 6 5];

highdoses = highdoses(resort);

figure;

plot(scores_early(2,~ismember(celltypes,highdoses) & ~outliers),scores_early(3,~ismember(celltypes,highdoses) & ~outliers),'o','MarkerSize',1,'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor','none')

hold on

color_ind = 1;
colmap = [linspace(0,1,length(highdoses)+1)' ones(length(highdoses)+1,1) ones(length(highdoses)+1,1)*.9];
colmap = hsv2rgb(colmap(1:end-1,:));
legstr = cell(length(highdoses),1);
legh = [];
markers = {'o','s','v','d','^','>'};
for isite = highdoses([6 2 3 4 5 1])
    sprop = siteprop(isite);
    legstr{isite == highdoses} = sprop.lig_name(1:3);
    
    scores = scores_early(:,celltypes == isite);
    legh(isite == highdoses) = plot(scores(2,~outliers(celltypes == isite)),scores(3,~outliers(celltypes == isite)),markers{isite == highdoses},'MarkerSize',5,'MarkerFaceColor',colmap(isite == highdoses,:),'MarkerEdgeColor','none');
    plotEllipsis(scores(2,~outliers(celltypes == isite)),scores(3,~outliers(celltypes == isite)),colmap(isite == highdoses,:),.5);
end


xlim = [-.25 .35];
ylim = [-.05 .15];
set(gca,'XLim',xlim,'YLim',ylim)

h = legend(legh,legstr);
ch = get(h,'child');
for ileg = 1:length(ch)/3
    ilegch = (ileg-1)*3+2;
    set(ch(ilegch),'LineStyle','-','LineWidth',1,'Color',colmap(size(colmap,1)-ileg+1,:)); 
end

ylabel('score of transient harmonic')
xlabel('score of sustained harmonic')
