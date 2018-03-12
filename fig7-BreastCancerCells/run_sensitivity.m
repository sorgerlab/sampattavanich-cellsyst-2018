addpath util

% Number of learning runs to perform with addition of noise 
ns = 100;

% List of cell lines
cell_lines = {'184A1', 'MCF10A', 'SKBR3', 'HCC1806', ...
            'HS578T', 'MDA231', 'BT20', 'MCF7', 'T47D'};

% Read data for each cell line
for c=1:length(cell_lines)
    data{c} = read_data_table(cell_lines{c});
end

% Create table for learning results
columns = {'cell_line', 'algo', 'pAE', 'pEA', 'pEF', 'pAF'};
sens = cell2table(cell(0, 6), 'VariableNames', columns);

% Run learning for each cell line
for c=1:length(cell_lines)
    fprintf('Sensitivities for %s...\n', cell_lines{c});
    for i=1:ns
		% Get marginal log-likelihoods with multiple methods, adding
		% some noise each time
        [LDBN_Gauss, LDBN_BGe, LDBN_BDe] = learn_structure(data{c}, true);
		% Derive edge probabilities from marginal likelihoods
        pdbde = get_edge_probs(LDBN_BDe);
        pdbga = get_edge_probs(LDBN_Gauss);
        pdbge = get_edge_probs(LDBN_BGe);
		% Set entries in table
        sens = [sens; cell2table({cell_lines{c}, ...
            'DBN_BDe', pdbde.pAE, pdbde.pEA, pdbde.pEF, pdbde.pAF}, ...
            'VariableNames', columns)];
        sens = [sens; cell2table({cell_lines{c}, ...
            'DBN_Gauss', pdbga.pAE, pdbga.pEA, pdbga.pEF, pdbga.pAF}, ...
            'VariableNames', columns)];
        sens = [sens; cell2table({cell_lines{c}, ...
            'DBN_BGe', pdbge.pAE, pdbge.pEA, pdbge.pEF, pdbge.pAF}, ...
            'VariableNames', columns)];
    end
end

% Plot Figure 7D
plot_edge_probs(sens(strcmp(sens.algo, 'DBN_Gauss'),:), cell_lines)
title('DBN learning')

% Plot Figure S7E
plot_edge_probs(sens(strcmp(sens.algo, 'DBN_BGe'),:), cell_lines)
title('DBN learning (BGe score)')

% Plot Figure S7F
plot_edge_probs(sens(strcmp(sens.algo, 'DBN_BDe'),:), cell_lines)
title('DBN learning (BDe score)')

% Print percentages reported in Figure S7D
means = varfun(@mean, sens, 'inputvariables', ...
    {'pAE', 'pEA', 'pEF', 'pAF'}, ...
    'groupingvariables', {'cell_line', 'algo'});

discr_prob = @(x) (x > 0.5);

means_discr = varfun(discr_prob, means, 'inputvariables', ...
    {'mean_pAE', 'mean_pEA', 'mean_pEF', 'mean_pAF'}, ...
    'groupingvariables', {'cell_line', 'algo'});

means_discr_bde = means_discr{strcmp(means_discr.algo, 'DBN_BDe'), ...
    {'Fun_mean_pAE', 'Fun_mean_pEA', 'Fun_mean_pEF', 'Fun_mean_pAF'}};
means_discr_gauss = means_discr{strcmp(means_discr.algo, 'DBN_Gauss'), ...
    {'Fun_mean_pAE', 'Fun_mean_pEA', 'Fun_mean_pEF', 'Fun_mean_pAF'}};
means_discr_bge = means_discr{strcmp(means_discr.algo, 'DBN_BGe'), ...
    {'Fun_mean_pAE', 'Fun_mean_pEA', 'Fun_mean_pEF', 'Fun_mean_pAF'}};
fprintf('BGe score vs BDe score structure overlap: %.0f%%\n', ...
		100*mean(mean(means_discr_bde ==  means_discr_bge)));
fprintf('BGe score vs Gaussian score structure overlap: %.0f%%\n', ...
		100*mean(mean(means_discr_gauss==  means_discr_bge)));
fprintf('BDe score vs Gaussian score structure overlap: %.0f%%\n', ...
		100*mean(mean(means_discr_bde ==  means_discr_gauss)));