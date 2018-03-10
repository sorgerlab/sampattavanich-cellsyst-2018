function plot_exampleEKAREVvsF3aN400()
[parentdir,~,~]=fileparts(pwd);
ndfilename = '03052014-r1.nd';
ndpathname = fullfile(parentdir,'rawdata','dualsensors');

prefix = ndfilename(1:(end-3));
[notp stagePos stageName channelnames] = readndfile(ndpathname,ndfilename);
outputsignalNo = 1;
sequenceNo1 = 1; %'Nuclei-EKAREV'
sequenceNo2 = 2; %'Cytosol-EKAREV'
sequenceNo3 = 3; %'Nuclei-F3aN400'
sequenceNo4 = 4; %'Cytosol-F3aN400'

sites = [14 5 19]; % 14 = BTC no drug, 5 = 4HRS addition of MEKi 1uM, 19 4HRs addition of AKTi 1UM
mytitle = {'G1:BTC,no drug'; 'G2:BTC,MEKi 1\muM at 4HR';'G3:BTC, AKTi 1\muM at 4HR'};
cellList{1} = [157 59];
cellList{2} = [261 222];
cellList{3} = [114 136];

for s = 1:length(sites)
    
    site = sites(s);
    tokens   = regexp(stageName{site}, 'r(?<row>\d+)c(?<col>\d+)|r(?<row>\d+)_c(?<col>\d+)|R(?<row>\d+)C(?<col>\d+)|R(?<row>\d+)_C(?<col>\d+)','tokens');
    if ~isempty(tokens)
        row = str2num(tokens{1}{1});
        col = str2num(tokens{1}{2});
    else
        row = site;
        col = 1;
    end
    field = 1;
    
    H5filename = ['H5OUT_r' num2str(row) '_c' num2str(col) '.h5'];
    signal_name = ['/field' num2str(field)  '/outputsignal' num2str(outputsignalNo)];
    timestamp_name = ['/field' num2str(field) '/timestamp' num2str(outputsignalNo)];
    selectedcells_name = ['/field' num2str(field) '/selectedcells'];
    
    warning off;
    
    if exist(fullfile(ndpathname,H5filename),'file')
        fileattrib(fullfile(ndpathname,H5filename),'+w');
        fid = H5F.open(fullfile(ndpathname,H5filename),'H5F_ACC_RDWR','H5P_DEFAULT');
        if H5L.exists(fid,signal_name,'H5P_DEFAULT') && ...
                H5L.exists(fid,timestamp_name,'H5P_DEFAULT') && ...
                H5L.exists(fid,selectedcells_name,'H5P_DEFAULT')
            H5F.close(fid);
            signalinfo = h5info(fullfile(ndpathname,H5filename), signal_name);
            countind = [signalinfo.Dataspace.Size(1) signalinfo.Dataspace.Size(2) 1];
            signal1 = permute(h5read(fullfile(ndpathname,H5filename),signal_name,double([1 1 sequenceNo1]), countind),[2 1 3]);
            signal2 = permute(h5read(fullfile(ndpathname,H5filename),signal_name,double([1 1 sequenceNo2]), countind),[2 1 3]);
            signal3 = permute(h5read(fullfile(ndpathname,H5filename),signal_name,double([1 1 sequenceNo3]), countind),[2 1 3]);
            signal4 = permute(h5read(fullfile(ndpathname,H5filename),signal_name,double([1 1 sequenceNo4]), countind),[2 1 3]);
            timestamp = h5read(fullfile(ndpathname,H5filename),timestamp_name);
            
            selected_cells = h5read(fullfile(ndpathname,H5filename),selectedcells_name);%   randi(size(signal1,2),1,min(size(signal1,2),200));
            
            c_signal = signal2; %Cytosolic EKAREV
            c_signal2 = signal4./signal3 ; %Cytosolic F3aN400/ Nuclear F3aN400
            c_signal(c_signal==0) = NaN;
            c_signal2(c_signal2==0) = NaN;
            c_signal(isinf(c_signal)) = NaN;
            c_signal2(isinf(c_signal2)) = NaN;
            
            i=1;
            for scell = cellList{s}
                figure(i);
                cX = [];
                cY = [];
                
                
                tpCriteria = ~isnan(c_signal(:,scell)) &  ~isnan(c_signal2(:,scell));
                pX = polyfit(timestamp(tpCriteria),c_signal(tpCriteria,scell),10);
                fX = polyval(pX,timestamp(tpCriteria));
                pY = polyfit(timestamp(tpCriteria),c_signal2(tpCriteria,scell),10);
                fY = polyval(pY,timestamp(tpCriteria));
                
                xtemp = c_signal(~isnan(c_signal(:,scell)),scell)- fX;% EKAREV detrended with polyfit
                ytemp = c_signal2(~isnan(c_signal2(:,scell)),scell)- fY; %F3aN400 C/N detrended  with polyfit
                
                cX = xtemp/(max(xtemp(timestamp<timestamp(66)))-min(xtemp(timestamp<timestamp(66)))); % EKAREV scaled to baseline range
                cY = ytemp/(max(ytemp(timestamp<timestamp(66)))-min(ytemp(timestamp<timestamp(66)))); % F3aN400 C/N scaled to baseline range
                
                goodTP = abs(cX) < 4 & abs(cY) < 4 ;
                cX = cX(goodTP);
                cY = cY(goodTP);
                goodtimestamp = timestamp(~isnan(c_signal2(:,scell)));
                goodtimestamp = goodtimestamp(goodTP);
                
                subplot(3,1,s);
                plot(goodtimestamp- timestamp(22),cX,'r');hold on; %recenter time of EKAREV based on ligand stimulation time
                plot(goodtimestamp- timestamp(22),cY,'b'); %recenter time of F3aN400 C/N based on ligand stimulation time
                
                ylim([-1 1]);
                xlim([timestamp(1)-timestamp(22) 800]);
                set(gca,'XTick',[0 200 400 600 800]);
                i=i+1;
                title(mytitle{s});
                legend({'EKAREV';'F3aN400 C/N'});
            end
        end
    end
end

function y = movingAverage(x, w)
k = ones(1, w) / w;
y = conv(x, k, 'same');
   
function y = movingIQR(x, w)
y = zeros(size(x));
for i = 1:(length(x)-w+1)
   cwin = x(i:(i+w-1));
   y(round(i+w/2)) = 1/iqr(cwin);
end
firstNonzero = find(y~=0,1,'first');
lastNonzero = find(y~=0,1,'last');
y(1:(firstNonzero-1)) = y(firstNonzero);
y((lastNonzero+1):end) = y(lastNonzero);




function [notp,stagePos,stageName,waveName] = readndfile(sourcefolder,filename)
% Search for number of string matches per line.
notp=-1;
stagePos = [];
stageName = [];
waveName = [];


if exist(fullfile(sourcefolder,filename),'file')
    fid = fopen(fullfile(sourcefolder,filename));
    y = 0;
    tline = fgetl(fid);
    sind = 1;
    wind = 1;
    notp=0;
    while ischar(tline)
        
        % Find number of time points
        
        testInd = regexp(tline,'NTimePoints');
        num = length(testInd);
        if num > 0
            tp  = regexp(tline, '(?<="NTimePoints", )\d+', 'match');
            notp = str2num(tp{1});
        end
        
        
        % Find stage naming
        testInd = regexp(tline,'Stage\d+');
        num = length(testInd);
        if num > 0
            stage  = regexp(tline, '(?<=")\w+(?=",)', 'match');
            stagePos{sind,1} = stage{1};
            stagename  = regexp(tline, '(?<="Stage\d+", ").+(?=")', 'match');
            stageName{sind,1} = stagename{1};
            sind=sind+1;
        end
        
        % Find stage naming
        testInd = regexp(tline,'WaveName\d+');
        num = length(testInd);
        if num > 0
            wavename1  = regexp(tline, '(?<="WaveName\d+", ")\w+(?=_)', 'match');
            wavename2  = regexp(tline, '(?<="WaveName\d+", "\w+_)\w+(?=")', 'match');
            wavename3  = regexp(tline, '(?<="WaveName\d+", ")\w+(?=")', 'match');
            if ~isempty(wavename1) && ~isempty(wavename2)
                waveName{wind} = ['w' num2str(wind) wavename1{1} '-' wavename2{1}];
            else
                waveName{wind} = ['w' num2str(wind) wavename3{1}];
            end
            wind=wind+1;
        end
        
        tline = fgetl(fid);
    end
    fclose(fid);
end

