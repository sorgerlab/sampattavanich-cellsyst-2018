% Figure 2B: Frequency Density Spectrum
addpath('./Functions/')

thres_sorted = [0 .1 .25 .5 .75 .9 1];

[parentdir,~,~]=fileparts(pwd);
load(fullfile(parentdir,'rawdata','Workspaces','fourier_signals_corrected_cleaned_newBTC2'));
time_range = [200 1475];

labels = {'< 10%', '10-25%', '25-50%', '50-75%', '75-90%', '> 90%', 'Pink noise 1', 'Pink noise 2', 'Pink noise 3', ...
    'Sine 1 + white noise 1', 'Sine 2 + white noise 1', 'Sine 3 + white noise 1', 'Sine 3 + white noise 2', 'Sine 3 + white noise 3'};
colmap = nan(length(c_signal_single),3);
colmap2 = gray(14);
colmap(1:6,:) = colmap2(4:9,:);
colmap(12,:) = [255,105,180]/255;
colmap(13,:) = [219,112,147]/255;
colmap(14,:) = [199,21,133]/255;
colmap(7,:) = [176,224,230]/255;
colmap(8,:) = [135,206,250]/255;
colmap(9,:) = [30,144,255]/255;
colmap(10,:) = [70,130,180]/255;
colmap(11,:) = [0,0,205]/255;

figure

linewidth = 2;

xfac = 1;
yfac = 1;
fontsize = 16;

setFigure(gcf,xfac,yfac,fontsize)
hold on

[~,range_ind_min] = min(abs(timestamp - time_range(1)));
[~,range_ind_max] = min(abs(timestamp - time_range(2)));
range_ind = range_ind_min:range_ind_max;

Fs = 1./((timestamp(2)-timestamp(1))*60); % Sampling every 5 min
L = length(range_ind);

NFFT = 2^nextpow2(L);

legh = [];
legstr = {};

for iplot = [1:6 12:14 7:11]

    c_signal_fft = [];
    c_signal_tmp = c_signal_single{iplot};
    c_signal_tmp(isinf(c_signal_tmp)) = nan;
    for i = 1:size(c_signal_tmp,2)
        if sum(isnan(c_signal_tmp(:,i))) > length(c_signal_tmp(:,i))-2
            c_signal_tmp(:,i) = 0;
        end
        if sum(isnan(c_signal_tmp(:,i))) > 0
            c_signal_tmp(:,i) = interp1(timestamp(~isnan(c_signal_tmp(:,i))),c_signal_tmp(~isnan(c_signal_tmp(:,i)),i),timestamp);
            vec = ~isnan(c_signal_tmp(:,i))';
            rl = find(vec ~= [vec(2:end), vec(end)+1]);
            data =  vec(rl);
            rl(2:end) = rl(2:end) - rl(1:end-1);
            if ~data(1)
                c_signal_tmp(1:rl(1),i) = c_signal_tmp(rl(1)+1,i);
            end
            if ~data(end)
                c_signal_tmp(end-rl(end)+1:end,i) = c_signal_tmp(end-rl(end),i);
            end
        end
            
        Y = fft(triang(length(range_ind)).*c_signal_tmp(range_ind,i),NFFT)/L;

        f = Fs/2*linspace(0,1,NFFT/2+1);

        c_signal_fft = [c_signal_fft 2*abs(Y(1:NFFT/2+1)).^2];
    end
    
    mean_fft = nanmean(c_signal_fft,2);
    std_fft = nanstd(c_signal_fft,[],2);
    std_fft = std_fft./sqrt(size(c_signal_fft,2)-1); % Std of mean --> 1/sqrt(N-1)
    
    legh = [legh semilogx(f,10*log10(mean_fft),'Color',colmap(iplot,:),'LineWidth',linewidth)];
end


set(gca,'XLim',[5e-5 1.2e-3])
set(gca,'XScale','log');
set(gca,'YLim',[-95 -60])

plot([1/(300*60) 1/(300*60)],get(gca,'YLim'),'k:','LineWidth',linewidth) % equals 300min period
plot([1/(15*60) 1/(15*60)],get(gca,'YLim'),'k:','LineWidth',linewidth) % equals 15min period

xlabel('Frequency (Hz)')
ylabel('log10 |Y(f)|')
ylabel('|Y(f)|^2 [dB]')

legend(legh,labels)

xlim3 = get(gca,'XLim');
h = axes('Position',get(gca,'Position'));
set(h,'XAxisLocation','top','Color','None','YTick',[],'XLim',xlim3,'XTick',[1/(300*60) 1/(80*60) 1/(15*60)],'XTickLabel',{'300','80','15'},'XScale','log')
xlabel('Period [min]')
