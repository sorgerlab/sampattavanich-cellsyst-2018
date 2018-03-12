function data = loadcsv(filename,datapath)

switch(filename)
    case {'130722_SCfeat','130722_Pav','140418_SCfeat','140418_Pav',...
            '140330_SCfeat','140330_Pav','140415_SCfeat','140415_Pav'}
        headerSpec = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s';
        formatSpec = '%d,%d,%d,%d,%s,%s,%s,%g,%g,%s,%g,%g,%g,%g,%g,%g,%g,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    case {'140215_SCdyn','140215_SCdyn_rev1'}
        headerSpec = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s';
        formatSpec = '%d,%d,%d,%d,%s,%s,%s,%g,%g,%d,%d,%d,%g,%g,%g,%g,%g,%g,%g,%g';
    case {'140324_Pav_rev1'}
        headerSpec = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s';
        formatSpec = '%d,%d,%d,%d,%s,%s,%s,%g,%g,%g,%g,%g,%g,%g,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    otherwise
        error(sprintf('Unknown dataset %s.\n',filename))
end

% filename = '130722_SCfeat.csv';
% filename = '130722_Pav.csv';
% filename = '140418_SCfeat.csv';
% filename = '140418_Pav.csv';
% filename = '140330_SCfeat.csv';
% filename = '140330_Pav.csv';
% headerSpec = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s';
% formatSpec = '%d,%d,%d,%d,%s,%s,%s,%g,%g,%s,%g,%g,%g,%g,%g,%g,%g,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';

% filename = '140415_SCfeat.csv';
% filename = '140415_Pav.csv';
% headerSpec = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s';
% formatSpec = '%d,%d,%d,%d,%s,%s,%s,%g,%g,%g,%g,%g,%g,%g,%d,%g,%g,%g,%g,%g,%g,%d,%g,%g,%g,%g';

% filename = '140215_SCdyn.csv';
% headerSpec = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s';
% formatSpec = '%d,%d,%d,%d,%s,%s,%s,%g,%g,%d,%d,%d,%g,%g,%g,%g,%g,%g,%g,%g';

headerSpec = strrep(headerSpec,',',' ');
formatSpec = strrep(formatSpec,',',' ');
formatSpec = strrep(formatSpec,'%g','%f');
formatSpec = strrep(formatSpec,'%d','%f');

fid = fopen(fullfile(datapath, [filename '.csv']));
hdrs = textscan(fid,headerSpec,1, 'delimiter',',');
tmpdata = textscan(fid,formatSpec,'delimiter',',');
fclose(fid);

outCell = cell(size(tmpdata{1},1), length(hdrs));
for i = 1:length(hdrs)
    if isnumeric(tmpdata{i})
        outCell(:,i) = num2cell(tmpdata{i});
    else
%         outCell(:,i) = tmpdata{i};
        outCell(:,i) = num2cell(nan(size(tmpdata{i})));
    end
end
data = cell2mat(outCell);


% Use the following to reproduce figure_3b_fPCA.m (130722_SCfeat):
% sites_all = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:69]; % Without FGF
% celltype = data(ismember(data(:,2),sites_all),2)';
% scores_early = data(ismember(data(:,2),sites_all),19:21)';
% Now run figure_3b_fPCA.m without the load('...') part

% Use the following to reproduce figure_3c_pulsatile.m (130722_SCfeat):
% sites_all = [4:10 17:-1:11 37:-1:31 44:50 57:-1:51 64:69]; % Without FGF
% celltypes = data(ismember(data(:,2),sites_all),2)';
% features = {'Final score','nEdge','SNR','Amplitude','Pairwise','Peak duration','Peak distance'};
% early_ws.scores_early = data(ismember(data(:,2),sites_all),19:21)';
% early_ws.celltypes = celltypes;
% scores_puls = data(ismember(data(:,2),sites_all),[29 25 26 24 5 27 28]);
% late_ws.celltypes = celltypes;
% late_ws.scores_puls = scores_puls;
% late_ws.features = features;
% Now run figure_3c_pulsatile.m without the load('...') part

% Use the following to reproduce figure_4a_fPCApulsing.m (140415_SCfeat):
% celltype = data(:,2)';
% dists = data(:,26)';
% scores_all = data(:,16:18)';
% Now run figure_4a_fPCApulsing.m without the load('...') part

% Use the following to reproduce figure_5d_pulsing_vs_early.m (140418_SCfeat):
% sites_egfmeki = [6:-1:1 7:12 19:24 18:-1:13 30:-1:25 31:39 41:72];
% egfmeki.celltype = data(ismember(data(:,2),sites_egfmeki),2)';
% egfmeki.dists = data(ismember(data(:,2),sites_egfmeki),29)';
% egfmeki_fPCA.scores_all = data(ismember(data(:,2),sites_egfmeki),19:21)';
% Now run figure_5d_pulsing_vs_early.m without the load('...') part

% Use the following to reproduce figure_5d_pulsing_vs_early.m (140330_SCfeat):
% sites_igfakti_unsorted = 1:72;
% sites_igfakti = sites_igfakti_unsorted([60:-1:1 61:72]);
% igfakti.celltype = data(ismember(data(:,2),sites_igfakti),2)';
% igfakti.dists = data(ismember(data(:,2),sites_igfakti),29)';
% igfakti.scores_all = data(ismember(data(:,2),sites_igfakti),19:21)';
% Now run figure_5d_pulsing_vs_early.m without the load('...') part

% Use the following to reproduce dualSensor_signalProcessing.m (140215_SCdyn)
% sites = [7	8	9   12	11	10  25	26	27  30	29	28  43	44	45  48	47	46  61	62	63];
% timestamp = data(1:209,13)+120;
% celltype = reshape(data(:,2),209,7833);
% celltype = celltype(1,:);
% c_signal_foxo = reshape(data(:,18),209,7833);
% c_signal_ekarev = reshape(data(:,15),209,7833);
% dualSensor_signalProcessing_checkCSV

% Calculate 
