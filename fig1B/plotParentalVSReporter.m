function plotParentalVSReporter(FigNo,FigName,tps,rows,cols,reporter_Mean,reporter_CIL,reporter_CIR,parental_Mean,parental_CIL,parental_CIR,mycolor)


c_reporter_Mean_combined = [];
c_parental_Mean_combined = [];

for t = tps
    c_reporter_Mean= reshape(reporter_Mean{t}(rows{t},cols{t}),[size(rows{t},2)*size(cols{t},2),1]);
    c_reporter_CIL  = reshape(reporter_CIL{t}(rows{t},cols{t}),[size(rows{t},2)*size(cols{t},2),1]);
    c_reporter_CIR  = reshape(reporter_CIR{t}(rows{t},cols{t}),[size(rows{t},2)*size(cols{t},2),1]);
    c_parental_Mean= reshape(parental_Mean{t}(rows{t},cols{t}),[size(rows{t},2)*size(cols{t},2),1]);
    c_parental_CIL  = reshape(parental_CIL{t}(rows{t},cols{t}),[size(rows{t},2)*size(cols{t},2),1]);
    c_parental_CIR  = reshape(parental_CIR{t}(rows{t},cols{t}),[size(rows{t},2)*size(cols{t},2),1]);
    
    figure(FigNo);errorbarxy(c_parental_Mean,c_reporter_Mean,c_parental_CIL,c_parental_CIR,c_reporter_CIL, c_reporter_CIR,{'k.', mycolor{t}, mycolor{t}}); hold on;
    figure(FigNo+2);plot(c_parental_Mean,c_reporter_Mean,'o','MarkerEdgeColor',mycolor{t}); hold on;
    
    c_reporter_Mean_combined = [c_reporter_Mean_combined;c_reporter_Mean];
    c_parental_Mean_combined = [c_parental_Mean_combined;c_parental_Mean];
    %xlim([-0.975 0]);
    %ylim([-1.15 0.2]);
end
rho  = corr(c_reporter_Mean_combined, c_parental_Mean_combined);
figure(FigNo);set(gcf,'Name',[FigName ' ' num2str(rho,2)]);
figure(FigNo+2);set(gcf,'Name',[FigName ' ' num2str(rho,2)]);
hold off;