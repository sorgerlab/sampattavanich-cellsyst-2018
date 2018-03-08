clear all; close all;
[parentdir,~,~]=fileparts(pwd);
load(fullfile(parentdir,'rawdata','combineddata06012014.mat'))
%%
% Fig. 1B
% Plot edogenous Endogenous FoxO3[log10(C_norm/N_norm)]  against F3aN400[log10(C_norm/N_norm)]
% Comparison of signals from time 15, 45, 90, 135
close all;
rows{1} = 8:-1:1;
rows{2} = 8:-1:1;
rows{3} = 8:-1:1; 
rows{4} = 8:-1:1;
rows{5} = 8:-1:1; 
cols{1} = 12:-1:7;
cols{2} = 12:-1:7;
cols{3} = 12:-1:7;
cols{4} = 12:-1:7;
cols{5} = 12:-1:7;
[reporter_Mean_EGF, reporter_CIL_EGF, reporter_CIR_EGF] = genFoxO3(0,reporter,rows,cols,'log10CoverN_4pixel_norm');
[parental_Mean_EGF, parental_CIL_EGF, parental_CIR_EGF] = genFoxO3(0,parental,rows,cols,'log10CoverN_4pixel_norm');
plotParentalVSReporter(1,'log10CoverN_4pixel_norm',2:5,rows,cols,reporter_Mean_EGF, reporter_CIL_EGF, reporter_CIR_EGF,parental_Mean_EGF, parental_CIL_EGF, parental_CIR_EGF,{'k','r','y','g','b'});
figure(3),hold on;plot([0 -2;0 2],[-1,0;1 0],':');xlim([-1 0.2]),ylim([-0.7 0.3]);
xlabel('Endogenous [log_1_0(normC/normN)]');
ylabel('F3aN400 [log_1_0(normC/normN)]');
figure(1),hold on;plot([0 -2;0 2],[-1,0;1 0],':');xlim([-1 0.2]),ylim([-0.7 0.3]);
xlabel('Endogenous [log_1_0(normC/normN)]');
ylabel('F3aN400 [log_1_0(normC/normN)]');


% %%
% % Plot raw version log10(C/N) using 4-pixel cytosolic mask
% rows{1} = 8:-1:1; %AKTi dose in uM = 0,0.0156,0.031,0.063,0.125,0.25,0.5,1
% rows{2} = 8:-1:1; %AKTi dose in uM = 0,0.0156,0.031,0.063,0.125,0.25,0.5,1
% rows{3} = 8:-1:1; %AKTi dose in uM = 0,0.0156,0.031,0.063,0.125,0.25,0.5,1
% rows{4} = 8:-1:1; %AKTi dose in uM = 0,0.0156,0.031,0.063,0.125,0.25,0.5,1
% rows{5} = 8:-1:1; %AKTi dose in uM = 0,0.0156,0.031,0.063,0.125,0.25,0.5,1
% cols{1} = 12:-1:7;%EGF dose in ng/mL = 0.4,1.2,3.7,11,33,100
% cols{2} = 12:-1:7;%EGF dose in ng/mL = 0.4,1.2,3.7,11,33,100
% cols{3} = 12:-1:7;%EGF dose in ng/mL = 0.4,1.2,3.7,11,33,100
% cols{4} = 12:-1:7;%EGF dose in ng/mL = 0.4,1.2,3.7,11,33,100
% cols{5} = 12:-1:7;%EGF dose in ng/mL = 0.4,1.2,3.7,11,33,100
% [reporter_Mean_EGF, reporter_CIL_EGF, reporter_CIR_EGF] = genFoxO3(0,reporter,rows,cols,'log10CoverN_4pixel');
% [parental_Mean_EGF, parental_CIL_EGF, parental_CIR_EGF] = genFoxO3(0,parental,rows,cols,'log10CoverN_4pixel');
% plotParentalVSReporter(1,'log10CoverN_4pixel',2:5,rows,cols,reporter_Mean_EGF, reporter_CIL_EGF, reporter_CIR_EGF,parental_Mean_EGF, parental_CIL_EGF, parental_CIR_EGF,{'k','r','y','g','b'});
% 
% %%
% % Comparison of signals from time 15, 45, 90, 135
% % Plot normalized version log10(C_norm/N_norm) using extended cytosolic mask
% rows{1} = 8:-1:1; %AKTi dose in uM = 0,0.0156,0.031,0.063,0.125,0.25,0.5,1
% rows{2} = 8:-1:1; %AKTi dose in uM = 0,0.0156,0.031,0.063,0.125,0.25,0.5,1
% rows{3} = 8:-1:1; %AKTi dose in uM = 0,0.0156,0.031,0.063,0.125,0.25,0.5,1
% rows{4} = 8:-1:1; %AKTi dose in uM = 0,0.0156,0.031,0.063,0.125,0.25,0.5,1
% rows{5} = 8:-1:1; %AKTi dose in uM = 0,0.0156,0.031,0.063,0.125,0.25,0.5,1
% cols{1} = 12:-1:7;%EGF dose in ng/mL = 0.4,1.2,3.7,11,33,100
% cols{2} = 12:-1:7;%EGF dose in ng/mL = 0.4,1.2,3.7,11,33,100
% cols{3} = 12:-1:7;%EGF dose in ng/mL = 0.4,1.2,3.7,11,33,100
% cols{4} = 12:-1:7;%EGF dose in ng/mL = 0.4,1.2,3.7,11,33,100
% cols{5} = 12:-1:7;%EGF dose in ng/mL = 0.4,1.2,3.7,11,33,100
% [reporter_Mean_EGF, reporter_CIL_EGF, reporter_CIR_EGF] = genFoxO3(0,reporter,rows,cols,'log10CoverN_extended_norm');
% [parental_Mean_EGF, parental_CIL_EGF, parental_CIR_EGF] = genFoxO3(0,parental,rows,cols,'log10CoverN_extended_norm');
% plotParentalVSReporter(2,'log10CoverN_extended_norm',2:5,rows,cols,reporter_Mean_EGF, reporter_CIL_EGF, reporter_CIR_EGF,parental_Mean_EGF, parental_CIL_EGF, parental_CIR_EGF,{'k','r','y','g','b'});
% %%
% % Comparison of signals from time 15, 45, 90, 135
% % Plot raw version using log10(C/N) using extended cytosolic mask
% rows{1} = 8:-1:1; %AKTi dose in uM = 0,0.0156,0.031,0.063,0.125,0.25,0.5,1
% rows{2} = 8:-1:1; %AKTi dose in uM = 0,0.0156,0.031,0.063,0.125,0.25,0.5,1
% rows{3} = 8:-1:1; %AKTi dose in uM = 0,0.0156,0.031,0.063,0.125,0.25,0.5,1
% rows{4} = 8:-1:1; %AKTi dose in uM = 0,0.0156,0.031,0.063,0.125,0.25,0.5,1
% rows{5} = 8:-1:1; %AKTi dose in uM = 0,0.0156,0.031,0.063,0.125,0.25,0.5,1
% cols{1} = 12:-1:7;%EGF dose in ng/mL = 0.4,1.2,3.7,11,33,100
% cols{2} = 12:-1:7;%EGF dose in ng/mL = 0.4,1.2,3.7,11,33,100
% cols{3} = 12:-1:7;%EGF dose in ng/mL = 0.4,1.2,3.7,11,33,100
% cols{4} = 12:-1:7;%EGF dose in ng/mL = 0.4,1.2,3.7,11,33,100
% cols{5} = 12:-1:7;%EGF dose in ng/mL = 0.4,1.2,3.7,11,33,100
% [reporter_Mean_EGF, reporter_CIL_EGF, reporter_CIR_EGF] = genFoxO3(0,reporter,rows,cols,'log10CoverN_extended');
% [parental_Mean_EGF, parental_CIL_EGF, parental_CIR_EGF] = genFoxO3(0,parental,rows,cols,'log10CoverN_extended');
% plotParentalVSReporter(5,'log10CoverN_extended',2:5,rows,cols,reporter_Mean_EGF, reporter_CIL_EGF, reporter_CIR_EGF,parental_Mean_EGF, parental_CIL_EGF, parental_CIR_EGF,{'k','r','y','g','b'});
% 
