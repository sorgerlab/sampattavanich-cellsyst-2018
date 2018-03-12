function [discenter,disL,disR] = genFoxO3(FigNo,mydataMatrix,Row,Col,log10CoverNtype)
%myxlim = [-1.6 0.6];
tp = fieldnames(mydataMatrix);
RSize = 8;
CSize = 6;
for t = 1:length(tp)
    i=1;
    for r = Row{t}
        for c = Col{t}
            c_Cytos = [];
            c_Nucs = [];

            FoxO3CpNS = mydataMatrix.(tp{t}){r,c}.(log10CoverNtype);
            currentMedian = prctile(FoxO3CpNS,50);
            currentP25 = prctile(FoxO3CpNS,25);
            currentP75 = prctile(FoxO3CpNS,75);
            
            discenter{t}(r,c) = currentMedian;%mean(FoxO3CpNS);
            myCI = paramci(fitdist(FoxO3CpNS,'Normal'),'Alpha',.01);
            disL{t}(r,c) = currentMedian-currentP25;%discenter{t}(r,c)-myCI(1,1);
            disR{t}(r,c) = currentP75-currentMedian; %myCI(2,1)-discenter{t}(r,c);
            % Plot histogram
            if FigNo>0
                figure(FigNo+t);subplot(RSize,CSize,i);
                [y,x] = hist(FoxO3CpNS,128);plot(x,y/max(y),'r-');hold on;
                plot([currentMedian currentMedian],[0 1],'k-'); hold on;
                plot([currentP25 currentP75],[0.5 0.5],'k-');
                title([num2str(r) ',' num2str(c)]);
                xlim([-1 1]);
                hold off;
                drawnow;
            end
            
            i=i+1;
        end
    end
end