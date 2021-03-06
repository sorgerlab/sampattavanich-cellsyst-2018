function comparedataat15Min()
[parentdir,~,~]=fileparts(pwd);
load(fullfile(parentdir,'rawdata','parentalVSReporter','combineddata06012014.mat'))

stepsize = 150;
% 100 ng/mL EGF
tp = fieldnames(reporter);
% plot reporter 100 ng/mL EGF
figure(1);subplot(2,1,1);
[y,x] = hist(reporter.(tp{2}){8,7}.log10CoverN_4pixel_norm,stepsize);
plot(x,y/max(y),'-','Color','r');hold on;
clear x y;

% plot parental 100 ng/mL EGF
[y,x] = hist(parental.(tp{2}){8,7}.log10CoverN_4pixel_norm,stepsize);
plot(x,y/max(y),'--','Color','r');
clear x y;

% plot reporter 100 ng/mL EGF with preincubation of 1uM of AKTi

[y,x] = hist(reporter.(tp{2}){1,7}.log10CoverN_4pixel_norm,stepsize);
plot(x,y/max(y),'-','Color','b');
clear x y;

% plot F3aN400 100 ng/mL EGF with preincubation of 1uM of AKTi
[y,x] = hist(parental.(tp{2}){1,7}.log10CoverN_4pixel_norm,stepsize);
plot(x,y/max(y),'--','Color','b');hold on;
xlabel('log_1_0(normC/normN)');
plot([0,0],[0,1],':k');hold off;

xlim([-2 1]);


figure(1);subplot(2,1,2);
[y,x] = hist(reporter.(tp{2}){8,7}.log10CoverN_4pixel,stepsize);
plot(x,y/max(y),'-','Color','r');hold on;
clear x y;

% plot parental 100 ng/mL EGF
[y,x] = hist(parental.(tp{2}){8,7}.log10CoverN_4pixel,stepsize);
plot(x,y/max(y),'--','Color','r');
clear x y;

% plot reporter 100 ng/mL EGF with preincubation of 1uM of AKTi

[y,x] = hist(reporter.(tp{2}){1,7}.log10CoverN_4pixel,stepsize);
plot(x,y/max(y),'-','Color','b');
clear x y;

% plot F3aN400 100 ng/mL EGF with preincubation of 1uM of AKTi
[y,x] = hist(parental.(tp{2}){1,7}.log10CoverN_4pixel,stepsize);
plot(x,y/max(y),'--','Color','b');hold on;

legend({'F3aN400:EGF','Endogenous:EGF','F3aN400:EGF+AKTi','Endogenous:EGF+AKTi'});

plot([0,0],[0,1],':k');hold off;
xlabel('log_1_0(C/N)');
xlim([-2 1]);
set(gcf,'Name','4-Pixel');

%%        

stepsize = 200;
% plot reporter 100 ng/mL EGF,extended
figure(2);
[y,x] = hist(reporter.(tp{2}){8,7}.log10CoverN_extended_norm,stepsize);
plot(x,y/max(y),'-','Color','r');hold on;
clear x y;

% plot reporter 100 ng/mL EGF,4pixel
[y,x] = hist(reporter.(tp{2}){8,7}.log10CoverN_4pixel_norm,stepsize);
plot(x,y/max(y),'--','Color','r');
clear x y;

% plot reporter 100 ng/mL EGF with preincubation of 1uM of AKTi, extended

[y,x] = hist(reporter.(tp{2}){1,7}.log10CoverN_extended_norm,stepsize);
plot(x,y/max(y),'-','Color','b');
clear x y;

% plot repoter 100 ng/mL EGF with preincubation of 1uM of AKTi, 4pixel
[y,x] = hist(reporter.(tp{2}){1,7}.log10CoverN_4pixel_norm,stepsize);
plot(x,y/max(y),'--','Color','b');hold on;

xlabel('log_1_0(normC/normN)');
plot([0,0],[0,1],':k');

xlim([-1.3 1]);

legend({'F3aN400-Extended:EGF','F3aN400-4Pixel:EGF','F3aN400-Extended:EGF+AKTi','F3aN400-4Pixel:EGF+AKTi'});

plot([0,0],[0,1],':k');hold off;
xlabel('log_1_0(C/N)');
set(gcf,'Name','Extended');

