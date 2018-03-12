function calculate_correlation()
[parentdir,~,~]=fileparts(pwd);
ndfilename = '03052014-r1.nd';
ndpathname = fullfile(parentdir,'rawdata','dualsensors');

prefix = ndfilename(1:(end-3));
[notp stagePos stageName channelnames] = readndfile(ndpathname,ndfilename);
outputsignalNo = 1;
%FOXO3a
%sequenceNo1 = 3;
%sequenceNo2 = 4;

%EKAREV
sequenceNo1 = 1;
sequenceNo2 = 2;
sequenceNo3 = 3;
sequenceNo4 = 4;
NoColoredLine = 12;
windowSize = 20;
sites = [14 5 19]; % 14 = BTC no drug, 5 = 4HRS addition of MEKi 1uM, 19 4HRs addition of AKTi 1UM
Mytitle = {'G1:BTC,no drug'; 'G2:BTC,MEKi 1\muM at 4HR';'G3:BTC, AKTi 1\muM at 4HR'};
cellList{1} = [157 16 136 55 395 231 18 3 357 207 145 333 127 234 116 335 127 376 188 412 219 144 109 165 493 94 347 59 215 171 126 105 113 442];
cellList{2} = [78 239 112 136 62 63 224 178 261 122 222 210 50 278 7];
cellList{3} = [146 34 269 136 97 140 114 101 55 38 17 210 304 171];


plotInd = 1;
GTime = [];
GRHO = [];
GType = [];

for s = 1:length(sites)
    
    site = sites(s);
    saveInd = 1;
    RHO = [];
    TIME = [];
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
            
            c_signal = signal2;%./signal1(chosenSegment,selected_cells(scell)) ;
            
            c_signal2 = signal4./signal3 ;
            c_signal(c_signal==0) = NaN;
            c_signal2(c_signal2==0) = NaN;
            c_signal(isinf(c_signal)) = NaN;
            c_signal2(isinf(c_signal2)) = NaN;
            
            
            for scell = cellList{s}
                cX = [];
                cY = [];
                
                if   length(find(~isnan(c_signal(:,scell)))) > 0.8*length(c_signal(:,scell))  && ...
                        length(find(~isnan(c_signal2(:,scell)))) > 0.8*length(c_signal2(:,scell))
                    
                    tpCriteria = ~isnan(c_signal(:,scell)) &  ~isnan(c_signal2(:,scell));
                    pX = polyfit(timestamp(tpCriteria),c_signal(tpCriteria,scell),5);
                    fX = polyval(pX,timestamp(tpCriteria));
                    pY = polyfit(timestamp(tpCriteria),c_signal2(tpCriteria,scell),5);
                    fY = polyval(pY,timestamp(tpCriteria));
                    %figure(1);
                    %plot(timestamp(tpCriteria),c_signal(tpCriteria,scell),'r');
                    %hold on;plot(timestamp(~isnan(c_signal(:,scell))),fX,'r:')
                    %plot(timestamp(tpCriteria),c_signal2(tpCriteria,scell),'b');
                    %plot(timestamp(tpCriteria),fY,'b:'); hold off;
                    
                    xtemp = c_signal(~isnan(c_signal(:,scell)),scell)- fX;%
                    ytemp = c_signal2(~isnan(c_signal2(:,scell)),scell)- fY; %
                    %figure(2);
                    
                    cX = xtemp/(max(xtemp(timestamp<timestamp(66)))-min(xtemp(timestamp<timestamp(66))));
                    cY = ytemp/(max(ytemp(timestamp<timestamp(66)))-min(ytemp(timestamp<timestamp(66))));
                    
                    goodTP = abs(cX) < 4 & abs(cY) < 4 ;
                    cX = cX(goodTP);
                    cY = cY(goodTP);
                    goodtimestamp = timestamp(~isnan(c_signal2(:,scell)));
                    goodtimestamp = goodtimestamp(goodTP);
                    
                    for w = 1:(size(signal1,1)-windowSize+1)
                        testTP = goodtimestamp>=timestamp(w) & goodtimestamp<timestamp((w+windowSize-1)) ;
                        
                        if ~isempty(testTP) && ~isempty(cX(testTP)) && ~isempty(cY(testTP))
                            temp = corr(cX(testTP),cY(testTP),'type','Pearson','rows','complete');
                            if ~isnan(temp) && ~isinf(temp)
                                RHO(saveInd) = temp;
                                TIME(saveInd) = timestamp(w) - timestamp(22);%(timestamp(w)+timestamp((w+windowSize-1)))/2;
                            end
                        end
                        saveInd=saveInd+1;
                    end
                end
            end
        end
    end

    figure(4);subplot(3,1,plotInd);
    out = scatplot(TIME,RHO,'squares',0.15,100,5,1,5);
    xlim([min(TIME) 600]);
    title(Mytitle{s});
    figure(3);subplot(3,1,plotInd);
    contour(out.xi,out.yi,out.zif,10);
    xlim([min(TIME) 600]);
    title(Mytitle{s});
    figure(5);subplot(3,1,plotInd);
    imagesc(out.xi(1,:),out.yi(:,1),out.zif);
    xlim([min(TIME) 600]);
    set(gca,'YDir','normal','YLim',[0 1],'YTick',[0 0.5 1]);
    title(Mytitle{s});
    drawnow;
    
    Group = zeros(size(TIME));
    Group(TIME <235 & TIME > 0) = 1 ; % MIN
    Group(TIME > 250 & TIME < 485) = 2;%MIN

    GRHO = [GRHO;RHO(:)];
    GTime = [GTime;Group(:)];
    cType = (s)*ones(size(TIME));
    GType = [GType;cType(:)];
    plotInd=plotInd+1;
    
end



figure(6);

plotGroup = GTime == 1 | GTime == 2;
boxplot(GRHO(plotGroup),[GType(plotGroup) GTime(plotGroup) ],...
        'factorseparator',[1 2],'medianstyle','line','whisker',0,'notch','on');
ylim([-0.5 1]);
title('G1:BTC no drug,G2:BTC w/MEKi,G3: BTC w/AKTi T1:0-235MIN,T2:250-485MIN');

h=findobj(gca,'tag','Outliers');
set(h,'Visible','off');
ylabel('Rho');
xlabel('Treatment Group/Temporal Group');

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

