% Figure 5D: fPC2 vs Pulse score
egf_dose = [4 20]; % Choose EGF dose here
igf_dose = [20 100]; % Same for IGF

pc = 2;

addpath('./Functions/')

[parentdir,~,~]=fileparts(pwd);

egfmeki = load(fullfile(parentdir,'rawdata','Workspaces','dists_04182014'));
egfmeki_fPCA = load(fullfile(parentdir,'rawdata','Workspaces','scores_04182014'));
igfakti = load(fullfile(parentdir,'rawdata','Workspaces','c_signal_03302014'));

extension_egfmeki = '04-18-2014';
extension_igfakti = '03-30-2014';

sites_egfmeki_unsorted = [1:39 41:72];
sites_egfmeki = [6:-1:1 7:12 19:24 18:-1:13 30:-1:25 31:39 41:72];
sites_igfakti_unsorted = 1:72;
sites_igfakti = sites_igfakti_unsorted([60:-1:1 61:72]);

puls_thres = [.2 .5; .3 .3; nan nan];

meki_doses = [.1 .1/4 .1/4^2 .1/4^3 .1/4^4 0];
akti_doses = [.5 .1 .025 .00625 0];

egf_all_doses = [0 .16 .8 4 20 100];
igf_all_doses = [0 .8 4 20 100];

cell_names = {'MCF10A','184A1','HCC1806'};
lig_names = {'EGF','IGF','FGF','HRG','HGF','EPR','BTC','NS'};

egfmeki_pulsing = [];
igfakti_pulsing = [];
egfmeki_early = [];
igfakti_early = [];
egfmeki_celltype = [];
igfakti_celltype = [];
egfmeki_drugdose = [];
igfakti_drugdose = [];
egfmeki_ligind = [];
igfakti_ligind = [];
egfmeki_ligdose = [];
igfakti_ligdose = [];

for i = 1:length(sites_egfmeki)

    isite = sites_egfmeki(i);
    s = siteprop(isite,extension_egfmeki);

    icell = strmatch(s.celltype,cell_names,'exact');
    idrug = s.drug_dose;
    ilig = strmatch(s.lig_name,lig_names,'exact');
    idose = s.lig_dose;

    egfmeki_pulsing = [egfmeki_pulsing sum(egfmeki.dists(egfmeki.celltype == isite) > puls_thres(icell,1)) / sum(egfmeki.celltype == isite)];
    egfmeki_early = [egfmeki_early nanmean(egfmeki_fPCA.scores_all(pc,egfmeki.celltype == isite))];
    egfmeki_celltype = [egfmeki_celltype icell];
    egfmeki_drugdose = [egfmeki_drugdose idrug];
    egfmeki_ligind = [egfmeki_ligind ilig];
    egfmeki_ligdose = [egfmeki_ligdose idose];

end
    
for i = 1:length(sites_igfakti)

    isite = sites_igfakti(i);
    s = siteprop(isite,extension_igfakti);

    icell = strmatch(s.celltype,cell_names,'exact');
    idrug = s.drug_dose;
    ilig = strmatch(s.lig_name,lig_names,'exact');
    idose = s.lig_dose;

    igfakti_pulsing = [igfakti_pulsing sum(igfakti.dists(igfakti.celltype == isite) > puls_thres(icell,1)) / sum(igfakti.celltype == isite)];
    igfakti_early = [igfakti_early nanmean(igfakti.scores_all(pc,igfakti.celltype == isite))];
    igfakti_celltype = [igfakti_celltype icell];
    igfakti_drugdose = [igfakti_drugdose idrug];
    igfakti_ligind = [igfakti_ligind ilig];
    igfakti_ligdose = [igfakti_ligdose idose];

end


cols{1,1} = [linspace(0,0.1,length(meki_doses)); linspace(0,0.6,length(meki_doses)); linspace(.1,1,length(meki_doses))]';
cols{1,2} = [linspace(0,0,length(meki_doses)); linspace(0,0.1,length(meki_doses)); linspace(.1,1,length(meki_doses))]';
cols{2,1} = [linspace(.1,1,length(akti_doses));linspace(0,0.6,length(akti_doses)); linspace(0,0.1,length(akti_doses))]';
cols{2,2} = [linspace(.1,1,length(akti_doses)); linspace(0,0.1,length(akti_doses)); linspace(0,0.1,length(akti_doses))]';

for icell = 1:2
    figure
	hold on
    
    for i = 1:length(sites_egfmeki)
        isite = sites_egfmeki(i);
        if egfmeki_celltype(i) == icell
            if any(egfmeki_ligdose(i) == egf_dose) && any(egfmeki_drugdose(i) == meki_doses)
                col = cols{1,egfmeki_ligdose(i) == egf_dose}(egfmeki_drugdose(i) == meki_doses,:);
                doplot = 0;
                if egfmeki_ligind(i) == strmatch('EGF',lig_names,'exact')
                    if ismember(egfmeki_ligdose(i),egf_dose)
                        doplot = 1;
                    end
                end
                if doplot

                    plot(egfmeki_early(i),egfmeki_pulsing(i),'d','MarkerSize',12,'MarkerFaceColor',col,'MarkerEdgeColor',cols{1,egfmeki_ligdose(i) == egf_dose}(end,:))
                end
            end
        end
    end
    
    for i = 1:length(sites_igfakti)
        isite = sites_igfakti(i);
        if igfakti_celltype(i) == icell
            if any(igfakti_ligdose(i) == igf_dose) && any(igfakti_drugdose(i) == akti_doses)
                col = cols{2,igfakti_ligdose(i) == igf_dose}(igfakti_drugdose(i) == akti_doses,:);
                doplot = 0;
                if igfakti_ligind(i) == strmatch('IGF',lig_names,'exact')
                    if ismember(igfakti_ligdose(i),igf_dose)
                        doplot = 1;
                    end
                end
                if doplot

                    plot(igfakti_early(i),igfakti_pulsing(i),'o','MarkerSize',12,'MarkerFaceColor',col,'MarkerEdgeColor',cols{2,igfakti_ligdose(i) == igf_dose}(end,:))
                end
            end
        end
    end

    title(cell_names{icell})
    set(gca,'XLim',[-.15 .25],'YLim',[0 1])
    xlabel(sprintf('Early PC%i',pc))
    ylabel('Fraction of pulsing cells')
    
end