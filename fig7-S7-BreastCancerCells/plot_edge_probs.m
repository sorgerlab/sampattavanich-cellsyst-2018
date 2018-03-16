function plot_edge_probs(sens, cellLines)
	% Plots the probabilities of edges between ERK, AKT and FOXO3a given
	% the results of DBN learning with a given scoring method,
	% and a list of cell lines
    figure; hold on;
    colors = {'k','r',[0.2,0.2,1],[0,0.7,0]};
    markers = {'v','^','o','d'};
    linewidth = 1.5;
    markerfacecolor = 'auto';
    markersize = 12;
    for c=1:length(cellLines)
        sens_cell_line = sens(strcmp(sens.cell_line,cellLines{c}),:);
        yloc = length(cellLines)-c+1;

        pAE = mean(sens_cell_line.pAE);
        pEA = mean(sens_cell_line.pEA);
        pEF = mean(sens_cell_line.pEF);
        pAF = mean(sens_cell_line.pAF);

        sAE = std(sens_cell_line.pAE);
        sEA = std(sens_cell_line.pEA);
        sEF = std(sens_cell_line.pEF);
        sAF = std(sens_cell_line.pAF);

        % Marker is at no noise point
        plot([pAE-sAE,pAE+sAE],[yloc,yloc]+0.1,'color',colors{1},'linewidth',2);
        plot([pEA-sEA,pEA+sEA],[yloc,yloc]+0.05,'color',colors{2},'linewidth',2);
        plot([pEF-sEF,pEF+sEF],[yloc,yloc]-0.05,'color',colors{3},'linewidth',2);
        plot([pAF-sAF,pAF+sAF],[yloc,yloc]-0.1,'color',colors{4},'linewidth',2);

        plot([0,1],[yloc,yloc],'k--');
        plot(pAE,yloc+0.1,markers{1},...
                'color',colors{1},...
                'markerfacecolor',markerfacecolor,...
                'markersize',markersize,...
                'linewidth',linewidth);
        plot(pEA,yloc+0.05,markers{2},...
                'color',colors{2},...
                'markerfacecolor',markerfacecolor,...
                'markersize',markersize,...
                'linewidth',linewidth);
        plot(pAF,yloc-0.1,markers{4},...
                'color',colors{4},...
                'markerfacecolor',markerfacecolor,...
                'markersize',markersize,...
                'linewidth',linewidth);
        plot(pEF,yloc-0.05,markers{3},...
                'color',colors{3},...
                'markerfacecolor',markerfacecolor,...
                'markersize',markersize,...
                'linewidth',linewidth);


    end
    xlim([0,1]);
    set(gca,'xtick',[0,0.5,1]);
    xlabel('Probability of edge');
    ylim([0,length(cellLines)+1]);
    set(gca,'ytick',1:length(cellLines));
    set(gca,'yticklabel',fliplr(cellLines));
    legends = {'AKT\rightarrow ERK',...
        'ERK\rightarrow AKT',...
        'ERK\rightarrow FoxO3a',...
        'AKT\rightarrow FoxO3a'};
end