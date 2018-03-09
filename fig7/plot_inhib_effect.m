function d_matrix = plot_inhib_effect()
	% Plot the effect of AKT inhibition on ERK and of ERK inhibition
	% on AKT.
    addpath util
    cell_lines = {'184A1','MCF10A','SKBR3','HCC1806','HS578T', ...
                'MDA231','BT20','MCF7','T47D'};
    for c = 1:length(cell_lines)
        data = read_data_table(cell_lines{c});
        ligands = unique(data.ligand);

        % ERK upon AKT inhibition
        akt_inhib = data((data.meki==0) & (data.akti==1),:);
        akt_noinhib = data((data.meki==0) & (data.akti==0),:);

        ts = akt_inhib.time(strcmp(akt_inhib.ligand, 'NS'));
        auct = @(x) auc(ts, x);

        akt_inhib_auc = varfun(auct, akt_inhib, ...
            'InputVariables', 'erk_mu', ...
            'GroupingVariables', 'ligand');

        akt_noinhib_auc = varfun(auct, akt_noinhib, ...
            'InputVariables', 'erk_mu', ...
            'GroupingVariables', 'ligand');

        akt_inhib_auc.diff = akt_inhib_auc.Fun_erk_mu - akt_noinhib_auc.Fun_erk_mu;


        % AKT upon MEK inhibition
        erk_inhib = data((data.akti==0) & (data.meki==1),:);
        erk_noinhib = data((data.akti==0) & (data.meki==0),:);

        erk_inhib_auc = varfun(auct, erk_inhib, ...
            'InputVariables', 'akt_mu', ...
            'GroupingVariables', 'ligand');

        erk_noinhib_auc = varfun(auct, erk_noinhib, ...
            'InputVariables', 'akt_mu', ...
            'GroupingVariables', 'ligand');

        erk_inhib_auc.diff = erk_inhib_auc.Fun_akt_mu - erk_noinhib_auc.Fun_akt_mu;


        % FOXO upon AKT inhibition
        foxo_akt_inhib = data((data.meki==0) & (data.akti==1),:);
        foxo_akt_noinhib = data((data.meki==0) & (data.akti==0),:);


        foxo_akt_inhib_auc = varfun(auct, foxo_akt_inhib, ...
            'InputVariables', 'foxo_pd', ...
            'GroupingVariables', 'ligand');

        foxo_akt_noinhib_auc = varfun(auct, foxo_akt_noinhib, ...
            'InputVariables', 'foxo_pd', ...
            'GroupingVariables', 'ligand');

        foxo_akt_inhib_auc.diff = foxo_akt_inhib_auc.Fun_foxo_pd - foxo_akt_noinhib_auc.Fun_foxo_pd;


        % FOXO upon ERK inhibition
        foxo_erk_inhib = data((data.meki==1) & (data.akti==0),:);
        foxo_erk_noinhib = data((data.meki==0) & (data.akti==0),:);

        foxo_erk_inhib_auc = varfun(auct, foxo_erk_inhib, ...
            'InputVariables', 'foxo_pd', ...
            'GroupingVariables', 'ligand');

        foxo_erk_noinhib_auc = varfun(auct, foxo_erk_noinhib, ...
            'InputVariables', 'foxo_pd', ...
            'GroupingVariables', 'ligand');

        foxo_erk_inhib_auc.diff = foxo_erk_inhib_auc.Fun_foxo_pd - foxo_erk_noinhib_auc.Fun_foxo_pd;


        for i = 1:length(ligands)
            dE(c,i) = akt_inhib_auc{strcmp(ligands{i}, akt_inhib_auc.ligand), {'diff'}};
            dA(c,i) = erk_inhib_auc{strcmp(ligands{i}, erk_inhib_auc.ligand), {'diff'}};
            dFA(c,i) = foxo_akt_inhib_auc{strcmp(ligands{i}, foxo_akt_inhib_auc.ligand), {'diff'}};
            dFE(c,i) = foxo_erk_inhib_auc{strcmp(ligands{i}, foxo_akt_inhib_auc.ligand), {'diff'}};
        end
    end

    clim = 60;

	% Plot Figure S7A and B, heat map plot of changes
    figure;
    colormap(get_colormap());
    subplot(1,2,1);
    imagesc(dA, [-clim, clim]);
    set(gca,'xticklabelrotation', 90)
    title('MEKi driven AKT change')
    colorbar;
    set(gca, 'yticklabel', cell_lines);
    set(gca, 'xticklabel', ligands);
    subplot(1,2,2);
    imagesc(dE, [-clim, clim]);
    set(gca,'xticklabelrotation', 90)
    title('AKTi-driven ERK change');
    colorbar;
    set(gca, 'yticklabel', cell_lines);
    set(gca, 'xticklabel', ligands);


	% Plot Figure 7C, same changes plotted with + markers
    figure;
    hold on;
	ylim([0 length(cell_lines)+1]);
    xlim([-clim clim]);
	set(gca,'Ydir','reverse')
	for i=1:length(cell_lines)
		p1 = plot(dE(i,:), i*ones(1, length(dE(i,:))), 'k+');
		p2 = plot(dA(i,:), i*ones(1, length(dA(i,:))), 'r+');
	end
	plot([0,0], [0,length(cell_lines)+1], 'k-');
	legend([p1 p2], {'\DeltaERK (w/ AKTi)', '\DeltaAKT (w/ MEKi)'}, ...
		   'location', 'southoutside');
    set(gca,'ytick',1:length(cell_lines));
    set(gca,'yticklabel',cell_lines);
	box on;
end

function cm = get_colormap()
    g = [linspace(0.3,1,32)';flipud(linspace(0.3,1,32)')];
    r = [linspace(0.3,1,32)'; ones(32,1)];
    b = [ones(32,1); flipud(linspace(0.3,1,32)')];
    cm = [r, g, b];
end
