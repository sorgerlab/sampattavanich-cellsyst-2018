% Figure 6C: ERK/FoxO3 pulse correlation
addpath('./functions/')

[parentdir,~,~]=fileparts(pwd);
datapath = fullfile(parentdir,'rawdata','Workspaces');
data = loadcsv('140215_SCdyn_rev1',datapath);

sites_all = [9 10 27 28 45 46 63];
extension = '02-15-2014_retracked';

timestamp = data(1:209,13)+120;
celltype = reshape(data(:,2),209,size(data,1)/209);
celltype = celltype(1,:);
c_signal_ekarev = reshape(data(:,19),209,size(data,1)/209);
c_signal_foxo = reshape(data(:,20),209,size(data,1)/209);

myCorrC = nan(1,size(c_signal_foxo,2));
for i = 1:length(myCorrC)
    corrInd = ~isnan(c_signal_foxo(:,i)) & ~isnan(c_signal_ekarev(:,i)) & timestamp > 200;
    if sum(corrInd) > 0
        mycor = corrcoef(c_signal_foxo(corrInd,i),c_signal_ekarev(corrInd,i));
        myCorrC(i) = mycor(1,2);
    end
end
[tmp,ind_corrC] = sort(myCorrC,'descend');

% Boxplot for high doses
highdoses = [9 10 27 28 45 46 63];
myext = '02-15-2014_retracked';

mylab = {};
boxCorrC = [];

for i = 1:length(highdoses)
    s = siteprop(highdoses(i),myext);
    
    mylab{end+1} = s.lig_name;
    boxCorrC = padconcatenation(boxCorrC,myCorrC(celltype == highdoses(i)),1);
    
end

figure
boxplot(boxCorrC')

set(gca,'XTick',1:length(highdoses),'XTickLabel',mylab)