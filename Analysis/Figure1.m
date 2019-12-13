clear
load('C:\Users\BuonoLab\Box Sync\Juan_2018\PV_SST_Opto_Study_19\ExcitabilityData\PyrPV\SingleExcitPYRPVData.mat')
Single = META;
clearvars -EXCEPT Single

load('C:\Users\BuonoLab\Box Sync\Juan_2018\PV_SST_Opto_Study_19\ExcitabilityData\PyrPV\PairedExcitPYRPVData.mat')
Paired = GOODPairedMETA;
clearvars -EXCEPT Single Paired

samplePV     = 10;
sampleWT     = 12;
AdaptCurrent = 8;

%% Getting Data from Individually Recorded Cells
% Getting Indexes of ProperVariables
T1 = find([Single.Group]==1); %PV Indexes
T2 = find([Single.Group]==2); %WT Indexes
% TargetInt1 = find(Single(T1(1)).CONDINT==[NaN NaN 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4]);
% TargetInt2 = find(Single(T2(1)).CONDINT==[NaN NaN 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4]);

%Calculating Variables for Spike Num
SingleSpikesPV = cell2mat({Single(T1).Spikenum}); SingleSpikesPV = SingleSpikesPV(2:10,:)';
SingleSpikesPYR = cell2mat({Single(T2).Spikenum}); SingleSpikesPYR = SingleSpikesPYR(2:10,:)';

%Calculating ISI for the example trace
ISIWT = diff(Single(sampleWT).SpData{8}(:,2)'); %from the sample trace at intensity 0.3nA
ISIPV = diff(Single(samplePV).SpData{8}(:,2)'); %from the sample trace at intensity 0.3nA

%Calculating Average Adaptation Index
AdaptIndexPV = cell2mat({Single(T1).AdaptIndex}); AdaptIndexPV = AdaptIndexPV(8,:)';
AdaptIndexPYR = cell2mat({Single(T2).AdaptIndex}); AdaptIndexPYR = AdaptIndexPYR(8,:)';
SingleAdaptIndex = {AdaptIndexPYR,AdaptIndexPV};

%% Getting Data from Simultaneous Recording
BadVarName    = squeeze(struct2cell(Paired)); 
Ind10Steps    = find(cellfun(@length,BadVarName(3,:))>9); %Finding Files with 10 or more Excitability current steps
Ind13Steps    = find(cellfun(@length,BadVarName(3,:))>12); %Finding Files with 13 or more Excitability current steps
Paired10Steps = Paired(Ind10Steps);
Paired13Steps = Paired(Ind13Steps);

PairedAdaptIndex      = zeros(length(Paired10Steps),2);
PairedSpikenumPYR     = zeros(length(Paired10Steps),9);
PairedSpikenumPV      = zeros(length(Paired13Steps),9);
PairedSpikenumPVLater = zeros(length(Paired13Steps),3);

for d = 1:length(Paired10Steps)
    PairedAdaptIndex(d,:)  = Paired10Steps(d).AdaptIndex(AdaptCurrent,:);
    PairedSpikenumPYR(d,:) = Paired10Steps(d).SpikenumPYR(2:10)'; %We don't count the first one because it is a negative current
    PairedSpikenumPV(d,:)  = Paired10Steps(d).SpikenumPV(2:10)';
end

for d = 1:length(Paired13Steps)
    PairedSpikenumPVLater(d,:) = Paired13Steps(d).SpikenumPV(11:13);
end

%% Combining Single and Simultaneous Data
AdaptIndex = {[SingleAdaptIndex{1};PairedAdaptIndex(:,1)],[SingleAdaptIndex{2};PairedAdaptIndex(:,2)]};

SpikesPV = [SingleSpikesPV;PairedSpikenumPV];
SpikesPYR = [SingleSpikesPYR;PairedSpikenumPYR];
SpikesPVLater = PairedSpikenumPVLater;

AvrSpikesPYR = mean(SpikesPYR);
SEMPYR = std(SpikesPYR)/sqrt(size(SpikesPYR,1));

%The N for the last 4 values in the SEMPV is different than the rest. 
AvrSpikesPV = [mean(SpikesPV) mean(SpikesPVLater)];
SEMPV =  [std(SpikesPV)/sqrt(size(SpikesPV,1)) std(SpikesPVLater)/sqrt(size(SpikesPVLater,1))];

XPYR = [0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45];
XPV = [XPYR 0.5 0.55 0.6];
%% Plotting Figure 1
close all

MainFigure = figure;
set(gcf,'color','w');
ScreenSize = get(groot,'ScreenSize');
MainFigure.Position = [0 0 ScreenSize(3)/1.1 ScreenSize(4)/1.1]; %This is to Force the figure to have a vertical orientation

SP1 = subplot('Position',[0.04 0.77 0.28 0.23]);
SP2 = subplot('Position',[0.36 0.77 0.28 0.23]);
SP3a = subplot('Position',[0.03 0.43 0.21 0.3 ]);
SP3b = subplot('Position',[0.2628 0.44 0.055 0.3]);
SP4a = subplot('Position',[0.38 0.43 0.21 0.3]);
SP4b = subplot('Position',[0.6128 0.44 0.055 0.3]);
SP5 = subplot('Position',[0.04 0.05 0.24 0.29]);
SP6 = subplot('Position',[0.38 0.05 0.24 0.29]);
AX1 = axes('Position',[0 0 1 1],'Visible','off');

SP7 = subplot('Position',[0.68 0.77 0.28 0.23]);
SP8a = subplot('Position',[0.72 0.43 0.21 0.3]);
SP8b = subplot('Position',[0.9428 0.44 0.055 0.3]);
SP9 = subplot('Position',[0.73 0.05 0.24 0.29]);
% AX2 = axes('Position',[0.09 0.43 0.296 0.3]);
% AX3 = axes('Position',[0.55 0.41 0.296 0.3]);

%%
% Plotting Raw Traces
subplot(SP1)
x1 = (0:3500)';
zWT = (Single(sampleWT).Data([1 3 5 8 9],1500:5000)');
xMat = repmat(x1,1,size(zWT,2));
zMatWT = zWT;
plot(xMat,zMatWT,'linewidth',2)
ylim([-70 50])
box off
axis off

%Crazy 3D stuff
% yMat = repmat((1:10),length(x1),1);
% plot3(xMat,yMat,zMat,'linewidth',2.5)
% view(40,40)
%print -djpeg100 PVTraces3D.jpg -r300

subplot(SP2)
zPV = (Single(samplePV).Data([1 3 5 8 9],1500:5000)');
zMatPV = zPV;
plot(xMat,zMatPV,'linewidth',2)
ylim([-70 30])
box off
axis off

%Plotting ISI and Adaptation Index
subplot(SP3a)
plot(1:length(ISIWT),ISIWT,'k-o','LineWidth',3,'MarkerFaceColor','k')
box off
set(gca,'fontsize',10,'fontweight','bold','XTick',[1 2 3 4 5])
title('Inter-Spike Interval for a PYR','FontSize',15)
SP3a.Title.Position = [3.5 73.7883 0];
xlabel('ISI Number','FontSize',12)
ylabel('ISI(ms)','FontSize',12)
xlim([0.75 length(ISIWT)+0.25])
ylim([0 70])

subplot(SP3b)
boxplot(AdaptIndex{1},'PlotStyle','compact','Colors','k','Labels',{''});
h = gca;
h.XAxis.Visible = 'off'; h.FontWeight = 'bold'; h.YAxisLocation = 'right';
% SP3b.Position = [0.36484375 0.43 0.05625 0.3];
% ylabel('Adaptation Idex','FontSize',12)
ylim([-0.025 0.15])
box off

subplot(SP4a)
plot(1:length(ISIPV),ISIPV,'r-o','LineWidth',3,'MarkerFaceColor','r')
box off
set(gca,'fontsize',10,'fontweight','bold')
title('Inter-Spike Interval for a PV','FontSize',15)
SP4a.Title.Position = [9.5 73.7883 0];
xlabel('ISI Number','FontSize',12)
ylabel('ISI(ms)','FontSize',12)
xlim([0.5 length(ISIPV)+0.25])
ylim([0 70])

subplot(SP4b)
boxplot(AdaptIndex{2},'PlotStyle','compact','Colors','r','Labels',{''})
h = gca;
% SP4b.Position = [0.84484375 0.43 0.05625 0.3];
h.XAxis.Visible = 'off'; h.FontWeight = 'bold'; h.YAxisLocation = 'right';
% ylabel('Adaptation Idex','FontSize',12)
ylim([-0.025 0.15])
box off

% Plotting SpikeNums
subplot(SP5)
errorbar(XPYR,AvrSpikesPYR,SEMPYR,'color','k','LineWidth',4,'Marker','o','MarkerSize',4)
set(gca,'fontsize',10,'fontweight','bold')
title('Input/Output PYR','FontSize',15)
xlabel('Input Intensity (nA)','FontSize',12)
ylabel('Spk Count','FontSize',12)
xlim([XPYR(1)*0.95 XPYR(end)*1.05])
ylim([0 18])
box off

subplot(SP6)
errorbar(XPV,AvrSpikesPV,SEMPV,'color','r','LineWidth',4,'Marker','o','MarkerSize',4)
hold on
% scatter([0.5 0.55 0.6],[14 18 20],55,[1 0 0],'d','MarkerFaceColor','r')
hold off
box off
set(gca,'fontsize',10,'fontweight','bold')
title('Input/Output PV','FontSize',15)
xlabel('Input Intensity (nA)','FontSize',12)
ylabel('Spk Count','FontSize',12)
xlim([XPV(1)*0.95 XPV(end)*1.05])
ylim([0 18])
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SST Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear AdaptIndex BadVarName Ind10Steps Ind13Steps  
load('C:\Users\BuonoLab\Box Sync\Juan_2018\PV_SST_Opto_Study_19\ExcitabilityData\PyrSST\PairedExcitSSTData.mat')
SSTMETA = GOODPairedMETA;

sampleSST = 2;
ISISST = diff(SSTMETA(sampleSST).spData{9,2}); %from the sample trace at intensity 0.3nA

BadVarName = squeeze(struct2cell(SSTMETA)); 
Ind10Steps = find(cellfun(@length,BadVarName(3,:))>9); %Finding Files with 10 or more Excitability current steps
% Ind13Steps = find(cellfun(@length,BadVarName(3,:))>12); %Finding Files with 13 or more Excitability current steps
SST10Steps = SSTMETA(Ind10Steps);
% SST13Steps = SSTMETA(Ind13Steps);

for j = 1:length(SST10Steps)
    AdaptIndex(j,:)   = SST10Steps(j).AdaptIndex(AdaptCurrent,2);
    SpikesSST(j,:)     = SST10Steps(j).SpikenumSST(2:10);
end

% for j = 1:length(SST13Steps) 
%     SpikesSSTLater(j,:) = SST13Steps(j).SpikenumSST(11:13);
% end

AvrSpikesSST = mean(SpikesSST);
SEMSST       =  std(SpikesSST)/sqrt(size(SpikesSST,1));
XSST         = [0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45];

%% Plotting SST

subplot(SP7)
zSST = (SSTMETA(sampleSST).Data{2}([1 3 5 8 9],6000:9500)');
zMatSST = zSST;
plot(xMat,zMatSST,'linewidth',2)
ylim([-76 36])
box off
axis off

subplot(SP8a)
plot(1:length(ISISST),ISISST,'b-o','LineWidth',3,'MarkerFaceColor','b')
box off
set(gca,'fontsize',10,'fontweight','bold')
title('Inter-Spike Interval for a SST','FontSize',15)
SP8a.Title.Position = [8.5 73.7883 0];
xlabel('ISI Number','FontSize',12)
ylabel('ISI(ms)','FontSize',12)
xlim([0.75 length(ISISST)+0.25])
ylim([0 70])

subplot(SP8b)
boxplot(AdaptIndex,'PlotStyle','compact','Colors','b','Labels',{''});
h = gca;
h.XAxis.Visible = 'off'; h.FontWeight = 'bold'; h.YAxisLocation = 'right';
% SP8b.Position = [0.36484375 0.43 0.05625 0.3];
% ylabel('Adaptation Idex','FontSize',12)
ylim([-0.025 0.15])
box off

subplot(SP9)
hold on
errorbar(XSST,AvrSpikesSST,SEMSST,'color','b','LineWidth',4,'Marker','o','MarkerSize',4)
% scatter([0.42 0.47 0.52],[2.97 4.35 6.1],55,[0 0 1],'d','MarkerFaceColor','b')
hold off
box off
set(gca,'fontsize',10,'fontweight','bold')
title('Input/Output SST','FontSize',15)
xlabel('Input Intensity (nA)','FontSize',12)
ylabel('Spk Count','FontSize',12)
xlim([XSST(1)*0.95 XSST(end)*1.05])
ylim([0 18])
box off

cd('C:\Users\BuonoLab\Desktop\Juan Luis\PV_SST_Opto_Study_19\FigureRepo')
print '-dtiffn' Figure1
