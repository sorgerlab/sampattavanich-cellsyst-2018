function [LDBN_Gauss,LDBN_BGe,LDBN_BDe] = learn_structure(data, add_noise)
    addpath ..
    addpath util

    if nargin < 2
        add_noise = false;
    end

    if add_noise
        noise_scale = 0.05;
        data.erk_mu = add_noise_to_signal(data.erk_mu, noise_scale);
        data.akt_mu = add_noise_to_signal(data.akt_mu, noise_scale);
        data.foxo_pd = add_noise_to_signal(data.foxo_pd, noise_scale);
    end

    % Prepare data for learning ERK incoming edges
    erk_table = data((data.meki==0), :);
    ligand_means = varfun(@mean, erk_table, 'inputvariables', 'erk_mu', ...
        'groupingvariables', 'ligand');
    erk_table.erk_mu_d = erk_table.erk_mu;
    erk_table.akt_mu_d = erk_table.akt_mu;
    erk_table.erk_mu_c = erk_table.erk_mu;
    for l=ligand_means.ligand'
        rows = strcmp(l, erk_table.ligand);
        erk_table.erk_mu_c(rows) = ...
            erk_table.erk_mu(rows) - ...
            ligand_means{strcmp(l, ligand_means.ligand), 'mean_erk_mu'};
        th = otsu(erk_table.erk_mu(rows));
        erk_table.erk_mu_d(rows) = ...
            1 * (erk_table.erk_mu(rows) > th);
    end
    th = otsu(erk_table.akt_mu);
    erk_table.akt_mu_d = 1 * (erk_table.akt_mu > th);

    Dminus = erk_table{erk_table.time ~= max(erk_table.time), ...
                {'erk_mu_c', 'akt_mu'}}';
    Dplus = erk_table{erk_table.time ~= min(erk_table.time), ...
                {'erk_mu_c', 'akt_mu'}}';
    Dminus_d = erk_table{erk_table.time ~= max(erk_table.time), ...
                {'erk_mu_d', 'akt_mu_d'}}';
    Dplus_d = erk_table{erk_table.time ~= min(erk_table.time), ...
                {'erk_mu_d', 'akt_mu_d'}}';
    V = ones(size(Dminus));
    S = 1;

    % Total nodes: AKT + ERK
    A = zeros(2);
    % ERK and AKT have 2 levels
    levels(1:2) = 2;

    for ae = [0,1]
        % Choose AKT to ERK edge
        A(2,1) = ae;
        lh_tmp = learn_dbn_bde(A,Dminus_d,Dplus_d,V,S,levels);
        LDBN_BDe.E(ae+1) = lh_tmp(1);
        lh_tmp = learn_dbn_gauss(A,Dminus,Dplus,V, false);
        LDBN_Gauss.E(ae+1) = lh_tmp(1);
        lh_tmp = learn_dbn_bge(A, Dminus, Dplus, V);
        LDBN_BGe.E(ae+1) = lh_tmp(1);
    end

    % Prepare data for learning AKT incoming edges
    akt_table = data(data.akti==0, :);
    ligand_means = varfun(@mean, akt_table, 'inputvariables', 'akt_mu', ...
        'groupingvariables', 'ligand');
    akt_table.akt_mu_d = akt_table.akt_mu;
    akt_table.akt_mu_c = akt_table.akt_mu;
    for l=ligand_means.ligand'
        akt_table.akt_mu_c(strcmp(l, akt_table.ligand)) = ...
            akt_table.akt_mu(strcmp(l, akt_table.ligand)) - ...
            ligand_means{strcmp(l, ligand_means.ligand), 'mean_akt_mu'};
        th = otsu(akt_table.akt_mu(strcmp(l, akt_table.ligand)));
        akt_table.akt_mu_d(strcmp(l, akt_table.ligand)) = ...
            1 * (akt_table.akt_mu(strcmp(l, akt_table.ligand)) > th);
    end
    th = otsu(akt_table.erk_mu);
    akt_table.erk_mu_d = 1 * (akt_table.erk_mu > th);

    Dminus = akt_table{akt_table.time ~= max(akt_table.time), ...
                {'erk_mu', 'akt_mu_c'}}';
    Dplus = akt_table{akt_table.time ~= min(akt_table.time), ...
                {'erk_mu', 'akt_mu_c'}}';
    Dminus_d = akt_table{akt_table.time ~= max(akt_table.time), ...
                {'erk_mu_d', 'akt_mu_d'}}';
    Dplus_d = akt_table{akt_table.time ~= min(akt_table.time), ...
                {'erk_mu_d', 'akt_mu_d'}}';

    V = ones(size(Dminus));
    S = 1;
    % Total nodes: AKT + ERK
    A = zeros(2);
    for ea = [0,1]
        % Choose AKT to ERK edge
        A(1,2) = ea;
        lh_tmp = learn_dbn_bde(A,Dminus_d,Dplus_d,V,S,levels);
        LDBN_BDe.A(ea+1) = lh_tmp(2);
        lh_tmp = learn_dbn_gauss(A,Dminus,Dplus,V, false);
        LDBN_Gauss.A(ea+1) = lh_tmp(2);
        lh_tmp = learn_dbn_bge(A, Dminus, Dplus, V);
        LDBN_BGe.A(ea+1) = lh_tmp(2);
    end

    % Prepare data for learning FOXO incoming edges
    foxo_table = data;
    % ERK, AKT and FOXO have 2 levels
    levels = [2, 2, 2, 2];
    th = otsu(foxo_table.erk_mu);
    foxo_table.erk_mu_d = 1 * (foxo_table.erk_mu > th);
    th = otsu(foxo_table.akt_mu);
    foxo_table.akt_mu_d = 1 * (foxo_table.akt_mu > th);

    th = otsu(foxo_table.foxo_pd);
    foxo_table.foxo_pd_d = 1 * (foxo_table.foxo_pd > th);

    Dminus = foxo_table{foxo_table.time ~= max(foxo_table.time), ...
                {'erk_mu', 'akt_mu', 'foxo_pd'}}';
    Dplus = foxo_table{foxo_table.time ~= min(foxo_table.time), ...
                {'erk_mu', 'akt_mu', 'foxo_pd'}}';
    Dminus_d = foxo_table{foxo_table.time ~= max(foxo_table.time), ...
                {'erk_mu_d', 'akt_mu_d', 'foxo_pd_d'}}';
    Dplus_d = foxo_table{foxo_table.time ~= min(foxo_table.time), ...
                {'erk_mu_d', 'akt_mu_d', 'foxo_pd_d'}}';

    V = ones(size(Dminus));
    S = 1;
    % Total nodes: AKT + ERK + FOXOmu + FOXOiqr
    A = zeros(3);
    for af = [0,1]
        for ef = [0,1]
            % Choose AKT to ERK edge
            A(1,3) = ef;
            A(2,3) = af;
            lh_tmp = learn_dbn_bde(A,Dminus_d,Dplus_d,V,S,levels);
            LDBN_BDe.F(2*af+ef+1) = lh_tmp(3);
            lh_tmp = learn_dbn_gauss(A,Dminus,Dplus,V, true);
            LDBN_Gauss.F(2*af+ef+1) = lh_tmp(3);
            lh_tmp = learn_dbn_bge(A, Dminus, Dplus, V);
            LDBN_BGe.F(2*af+ef+1) = lh_tmp(3);
        end
    end
end

function data = add_ligand_columns(data, type)
    ligands = unique(data.ligand);
    switch type
        case 'multiple'
            for i=1:length(ligands)
                if strcmp(ligands{i}, 'NS') == 1
                    continue
                else
                    data{:,ligands{i}} = 1*(strcmp(data.ligand,ligands{i}));
                end
            end
        case 'single'
            data{:, 'ligand_val'} = zeros(length(data.ligand), 1);
            for i=1:length(ligands)
                rows = (strcmp(data.ligand, ligands{i})==1);
                values = repmat(i-1, sum(rows), 1);
                data{rows, 'ligand_val'} = values;
            end
    end
end

function signal = add_noise_to_signal(signal, noise_scale)
    scale_factor = noise_scale * median(abs(signal(:)));
    signal = signal + scale_factor * randn(size(signal));
end
