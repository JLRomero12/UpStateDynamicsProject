clear

cd('C:\Users\BuonoLab\Desktop\Juan Luis\MainStudy\SOM-tdTomato')

load ExcitabilityMETA
Paired = META;
clear META

sampleCell = 2;

%% Getting Data
%Calculating ISI for the example trace
ISIPYR = diff(Paired(sampleCell).spData{9,1}); %from the sample trace at intensity 0.3nA
ISISOM = diff(Paired(sampleCell).spData{9,2}); %from the sample trace at intensity 0.3nA

%Getting Data from Structure into mat
AdaptIndex = zeros(length(Paired),2);
Spikes     = zeros(2,length(Paired),7);
SpikesLater= zeros(2,length(Paired)-2,3);
for j = 1:length(Paired)
    AdaptIndex(j,:)   = Paired(j).AdaptIndex(8,:);
    Spikes(1,j,:)     = Paired(j).SpikenumPYR(2:8,:)';
    Spikes(2,j,:)     = Paired(j).SpikenumPV(2:8,:)';
    if j>2
        SpikesLater(1,j,:)= Paired(j).SpikenumPYR(9:11,:)';
        SpikesLater(2,j,:)= Paired(j).SpikenumPV(9:11,:)';
    end
end

%% Combining Data and doing averages
%The N for the last 3 values in the SEMPV is different than the rest.
AvrSpikesPYR = [mean(squeeze(Spikes(1,:,:))) mean(squeeze(SpikesLater(1,:,:)))];
AvrSpikesSOM = [mean(squeeze(Spikes(2,:,:))) mean(squeeze(SpikesLater(2,:,:)))];
SEMPYR = [std(squeeze(Spikes(2,:,:)))/sqrt(size(squeeze(Spikes(2,:,:)),1)) std(squeeze(SpikesLater(2,:,:)))/sqrt(size(squeeze(SpikesLater(2,:,:)),1))];
SEMSOM =  [std(squeeze(Spikes(2,:,:)))/sqrt(size(squeeze(Spikes(2,:,:)),1)) std(squeeze(SpikesLater(2,:,:)))/sqrt(size(squeeze(SpikesLater(2,:,:)),1))];

% AdaptIndex = nanmean(AdaptIndex);

XPYR = [0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5];
XSOM = XPYR;
%% Plotting Figure 1
close all

MainFigure = figure;
set(gcf,'color','w');
ScreenSize = get(groot,'ScreenSize');
MainFigure.Position = [0 0 ScreenSize(3)/1.5 ScreenSize(4)/1.2]; %This is to Force the figure to have a vertical orientation

SP1 = subplot('Position',[0.03 0.77 0.43 0.23]);
SP2 = subplot('Position',[0.53 0.77 0.43 0.23]);
SP3a = subplot('Position',[0.08 0.43 0.28 0.3 ]);
SP3b = subplot('Position',[0.3648 0.37 0.05625 0.3]);
SP4a = subplot('Position',[0.56 0.43 0.28 0.3]);
SP4b = subplot('Position',[0.847 0.37 0.05625 0.3]);
SP5 = subplot('Position',[0.08 0.05 0.35 0.29]);
SP6 = subplot('Position',[0.56 0.05 0.35 0.29]);
AX1 = axes('Position',[0 0 1 1],'Visible','off');
% AX2 = axes('Position',[0.09 0.43 0.296 0.3]);
% AX3 = axes('Position',[0.55 0.41 0.296 0.3]);


% Plotting Raw Traces
subplot(SP1)
x1 = (0:3500)';
zWT = (Paired(sampleCell).Data{1}([1 3 5 8 9],1500:5000)');
xMat = repmat(x1,1,size(zWT,2));
zMatWT = zWT;
plot(xMat,zMatWT,'linewidth',2)
ylim([-90 20])
box off
axis off

%Crazy 3D stuff
% yMat = repmat((1:10),length(x1),1);
% plot3(xMat,yMat,zMat,'linewidth',2.5)
% view(40,40)
%print -djpeg100 PVTraces3D.jpg -r300

subplot(SP2)
zSOM = (Paired(sampleCell).Data{2}([1 3 5 8 9],6000:9500)');
zMatSOM = zSOM;
plot(xMat,zMatSOM,'linewidth',2)
ylim([-80 30])
box off
axis off

%Plotting ISI and Adaptation Index
subplot(SP3a)
plot(1:length(ISIPYR),ISIPYR,'k-o','LineWidth',3,'MarkerFaceColor','k')
box off
set(gca,'fontsize',10,'fontweight','bold')
title('Inter-Spike Interval for a PYR','FontSize',15)
SP3a.Title.Position = [3.5 70.7883 0];
xlabel('ISI Number','FontSize',12)
ylabel('ISI(ms)','FontSize',12)
xlim([0.75 length(ISIPYR)+0.25])
ylim([0 70])

subplot(SP3b)
boxplot(AdaptIndex(:,1),'PlotStyle','compact','Colors','k','Labels',{''});
h = gca;
h.XAxis.Visible = 'off'; h.FontWeight = 'bold'; h.YAxisLocation = 'right';
SP3b.Position = [0.36484375 0.43 0.05625 0.3];
ylabel('Adaptation Idex','FontSize',12)
ylim([-0.025 0.15])
box off

subplot(SP4a)
plot(1:length(ISISOM),ISISOM,'r-o','LineWidth',3,'MarkerFaceColor','r')
box off
set(gca,'fontsize',10,'fontweight','bold')
title('Inter-Spike Interval for a SOM','FontSize',15)
SP4a.Title.Position = [9.5 70.7883 0];
xlabel('ISI Number','FontSize',12)
ylabel('ISI(ms)','FontSize',12)
xlim([0.5 length(ISISOM)+0.25])
ylim([0 70])

subplot(SP4b)
boxplot(AdaptIndex(:,2),'PlotStyle','compact','Colors','r','Labels',{''})
h = gca;
SP4b.Position = [0.84484375 0.43 0.05625 0.3];
h.XAxis.Visible = 'off'; h.FontWeight = 'bold'; h.YAxisLocation = 'right';
ylabel('Adaptation Idex','FontSize',12)
ylim([-0.025 0.15])
box off

% Plotting SpikeNums
subplot(SP5)
errorbar(XPYR,AvrSpikesPYR,SEMPYR,'color','k','LineWidth',4,'Marker','o','MarkerSize',4)
set(gca,'fontsize',10,'fontweight','bold')
text([0.4 0.45 0.5],[15 18 20],'\ast','Color','k','FontWeight','bold','FontSize',16)
text(0.53,2,'\ast','Color','k','FontWeight','bold','FontSize',16)
text(0.55,2,['n=' num2str(size(SpikesLater(1,:,:),2)-2)],'Color','k','FontWeight','bold','FontSize',12)
title('Input/Output PYR','FontSize',15)
xlabel('Input Intensity (nA)','FontSize',12)
ylabel('Spk Count','FontSize',12)
xlim([XPYR(1)*0.95 XPYR(end)*1.05])
ylim([0 25])
box off

subplot(SP6)
errorbar(XSOM,AvrSpikesSOM,SEMSOM,'color','r','LineWidth',4,'Marker','o','MarkerSize',4)
box off
text([0.4 0.45 0.5],[15 18 20],'\ast','Color','r','FontWeight','bold','FontSize',16)
text(0.53,2,'\ast','Color','r','FontWeight','bold','FontSize',16)
text(0.55,2,['n=' num2str(size(SpikesLater(1,:,:),2)-2)],'Color','k','FontWeight','bold','FontSize',12)
set(gca,'fontsize',10,'fontweight','bold')
title('Input/Output SOM','FontSize',15)
xlabel('Input Intensity (nA)','FontSize',12)
ylabel('Spk Count','FontSize',12)
xlim([XSOM(1)*0.95 XSOM(end)*1.05])
ylim([0 25])
box off

axes(AX1)
text(0.38,0.42,['n=' num2str(length(AdaptIndex(~isnan(AdaptIndex(:,1)),1)))],'FontWeight','bold','FontSize',12)
text(0.855,0.42,['n=' num2str(length(AdaptIndex(~isnan(AdaptIndex(:,2)),1)))],'FontWeight','bold','FontSize',12)

[h,p,s,c] = ttest2(AdaptIndex{1}(~isnan(AdaptIndex{1})),AdaptIndex{2}(~isnan(AdaptIndex{2})));

% axes(AX2)
% plot(x1',Single(sampleWT).Data(8,1500:5000),'LineWidth',2,'Color',[0 0 0 0.3])
% box off
% axis off
%
% axes(AX3)
% plot(x1',Single(samplePV).Data(8,1500:5000),'LineWidth',2,'Color',[0 0 0 0.3])
% box off
% axis off