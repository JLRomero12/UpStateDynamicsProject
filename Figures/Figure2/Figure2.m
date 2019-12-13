clear
close all

ExamplePVCell = 'C2_0320_.mat';
UpStatePVTrace = [12 15 20];

ExampleSSTCell = 'S2C2_0518_.mat';
UpStateSSTTrace = [11];

GoodGreen = [0 0.7 0];
GoodCyan = [0 0.9 0.9];
%% INIT Graphics
MainFigure = figure;
% MainFigure = figure('WindowState','Maximized');
set(gcf,'color','w');
ScreenSize = get(groot,'ScreenSize');
MainFigure.Position = [0 0 ScreenSize(3)/1.15 ScreenSize(4)/1.4]; %This is to Force the figure to have a vertical orientation

% SP1  = subplot('Position',[0.05 0.85 0.925 0.15]);
% SP2  = subplot('Position',[0.05 0.7 0.925 0.15]);
% SP3  = subplot('Position',[0.05 0.53 0.925 0.15]);
% SP4  = subplot('Position',[0.05 0.38 0.925 0.15]);
% SP5  = subplot('Position',[0.05 0.05 0.45 0.3]);
% SP6  = subplot('Position',[0.6 0.05 0.35 0.3]);

SP1  = subplot('Position',[0.005 0.7505 0.68 0.22]);
SP2  = subplot('Position',[0.005 0.5005 0.68 0.22]);
SP3  = subplot('Position',[0.005 0.2505 0.68 0.22]);
SP4  = subplot('Position',[0.005 0.005 0.68 0.22]);
SP5  = subplot('Position',[0.75 0.52 0.23 0.45]);
SP6  = subplot('Position',[0.75 0.015 0.23 0.4]);


% SP11 = subplot('Position',[0.91 0.8750 0.08 0.125]);
% SP22 = subplot('Position',[0.91 0.725 0.08 0.125]);
% SP33 = subplot('Position',[0.91 0.535 0.08 0.125]);
% SP44 = subplot('Position',[0.91 0.385 0.08 0.125]);

%% Panel A   EXAMPLE PYR-PV

cd('..\..\ConnectivityExperiments\Pyr-PV-Recordings')
load(ExamplePVCell)

UpStateIndex = find(strcmp(CondOrder,'BMI_LIGHT'));
RawUpState1 = CH1(UpStateIndex);
RawUpState2 = CH2(UpStateIndex);
   
UpState1 = [];
UpState2 = [];

for i = 1:length(RawUpState1)
    upstate1 = RawUpState1{i};
    upstate2 = RawUpState2{i};
    UpState1 = [UpState1; upstate1];
    UpState2 = [UpState2; upstate2];    
end

dt = SampleRate(UpStateIndex(1));
x = (1:size(UpState1,2))*dt;

subplot(SP1)
plot(x,UpState1(UpStatePVTrace(2),:),'Color',GoodGreen,'LineWidth',2)
hold on
line([24000 24000],[-40 -10],'Color','k','LineWidth',2)
line([23000 24000],[-40 -40],'Color','k','LineWidth',2)
text(23125,-44,'1sec','FontSize',12,'FontWeight','bold')
h1 = text(24175,-14,'30mV','FontSize',12,'FontWeight','bold');
set(h1,'Rotation',270);
set(SP1,'xcolor','none')
set(SP1,'ycolor','none')
xlim([x(1) x(end)+1500])
ylim([-83 0])
[d1, X] = histcounts(UpState1(UpStatePVTrace(2),:));
plot(x(end)+400+d1,conv(X, [0.5 0.5], 'valid'),'Color',GoodGreen,'LineWidth',2)
% axis tight
axis off

% subplot(SP11)
% h = histogram([UpState1(UpStatePVTrace(1),:) UpState1(UpStatePVTrace(2),:) UpState1(UpStatePVTrace(3),:)]);
% hh = plot(conv(h.BinEdges, [0.5 0.5], 'valid'),h.BinCounts,'Color',GoodGreen,'LineWidth',1);
% set(gca,'xdir','reverse')
% camroll(270)
% xlim([-70 -35])
% ylim([0 10000])
% axis off

subplot(SP2)
plot(x,UpState2(UpStatePVTrace(2),:),'Color','r','LineWidth',2)
hold on
set(SP2,'xcolor','none')
set(SP2,'ycolor','none')
xlim([x(1) x(end)+1500])
ylim([-83 0])
[d2, X] = histcounts(UpState2(UpStatePVTrace(2),:));
plot(x(end)+400+d2,conv(X, [0.5 0.5], 'valid'),'Color','r','LineWidth',2)
% axis tight
axis off

% subplot(SP22)
% h = histogram([UpState2(UpStatePVTrace(1),:) UpState2(UpStatePVTrace(2),:) UpState2(UpStatePVTrace(3),:)]);
% hh = plot(conv(h.BinEdges, [0.5 0.5], 'valid'),h.BinCounts,'r','LineWidth',1);
% set(gca,'xdir','reverse')
% camroll(270)
% xlim([-70 -35])
% ylim([0 10000])
% axis off
%% Panel B EXAMPLE PYR-SST
%Possible Examples S1C1_0221 Condition 4, frame 1
%Possible Examples S1C1_0221 Condition 7, frame 5 6
%Possible Examples S1C1_0517 Condition 3, frame 4 5 6 Meh
%Possible Examples S1C2_0517 Condition 6, frame 2 3 5 6
%Possible Examples S2C2_0518 Condition 3, frame 4 5 6
%Possible Examples S3C1_0517 Condition 4, frame 5 6
%Possible Examples S3C1_0517 Condition 5, frame 2 3 5 6

cd('..\Pyr-SST-Recordings\')
load(ExampleSSTCell)

UpStateIndex = find(strcmp(CondOrder,'SPONT'));
RawUpState1 = CH1(UpStateIndex);
RawUpState2 = CH2(UpStateIndex);
   
UpState1 = [];
UpState2 = [];

for i = 1:length(RawUpState1)
    upstate1 = RawUpState1{i};
    upstate2 = RawUpState2{i};
    UpState1 = [UpState1; upstate1];
    UpState2 = [UpState2; upstate2];    
end

dt = SampleRate{UpStateIndex(1)};
x = (1:50000)*dt;

subplot(SP3)
plot(x,UpState1(12,55001:105000),'Color',GoodGreen,'LineWidth',2)
hold on
line([24000 24000],[-40 -10],'Color','k','LineWidth',2)
line([23000 24000],[-40 -40],'Color','k','LineWidth',2)
text(23125,-44,'1sec','FontSize',12,'FontWeight','bold')
h1 = text(24175,-16,'30mV','FontSize',12,'FontWeight','bold');
set(h1,'Rotation',270);
set(SP1,'xcolor','none')
set(SP1,'ycolor','none')
xlim([x(1) x(end)+1500])
ylim([-75 -5])
[d3, X] = histcounts(UpState1(12,55001:105000));
plot(x(end)+400+d3,conv(X, [0.5 0.5], 'valid'),'Color',GoodGreen,'LineWidth',2)
% axis tight
axis off
%8 is the previous one, possible are 13 and 16

% subplot(SP33)
% h = histcounts(UpState1(12,50001:110000));
% hh = plot(h.BinCounts,conv(h.BinEdges, [0.5 0.5], 'valid'),'Color',GoodGreen,'LineWidth',1);
% set(gca,'xdir','reverse')
% camroll(270)
% % xlim([-72 -37])
% % ylim([0 10000])
% axis off

subplot(SP4)
plot(x,UpState2(12,55001:105000),'Color',GoodCyan,'LineWidth',2)
hold on
set(SP2,'xcolor','none')
set(SP2,'ycolor','none')
xlim([x(1) x(end)+1500])
ylim([-75 -5])
[d4, X] = histcounts(UpState2(12,55001:105000));
plot(x(end)+400+d4,conv(X, [0.5 0.5], 'valid'),'Color',GoodCyan,'LineWidth',2)
% axis tight
axis off

% subplot(SP44)
% h = histogram([UpState2(UpStateSSTTrace(1),:) UpState2(UpStateSSTTrace(2),1:(end/4))]);
% hh = plot(conv(h.BinEdges, [0.5 0.5], 'valid'),h.BinCounts,GoodCyan,'LineWidth',1);
% set(gca,'xdir','reverse')
% camroll(270)
% xlim([-70 -35])
% ylim([0 10000])
% axis off
%% Panel C
cd('..\..\Figures\Figure2')

load('Pyr-ChR2PV-Paired-MasterMETA.mat')
PyrPVMeta = MasterMETA;
clear MasterMETA

load('SSTTdTomatoMETA.mat')
PyrSSTMeta = MasterMETA;
clear MasterMETA

load('1Trans-PVCheta-PairedMETA.mat')
OneTPyrPVMETA = MasterMETA;
clear MasterMETA
load('1Trans-PVHalo-PairedMETA.mat')
OneTPyrPVMETA = [OneTPyrPVMETA MasterMETA];
clear MasterMETA

load('1Trans-SSTCheta-PairedMETA.mat')
OneTPyrSSTMETA = MasterMETA;
clear MasterMETA
load('1Trans-SSTHalo-PairedMETA.mat')
OneTPyrSSTMETA = [OneTPyrSSTMETA MasterMETA];
clear MasterMETA

%Getting the Non-Stimulated Variables
PyrPVFR = nan(2,length(PyrPVMeta));
PyrPVFRMean = nan(2,length(PyrPVMeta));
for i = 1:length(PyrPVMeta)
    NonStimFlag = PyrPVMeta(i).StimFlag(PyrPVMeta(i).StimFlag ~= -1) == 0; % The FiringRate of Upstates labeled -1 was not calculated so if the StimFlag has -1 the dimensions won't match with the FR
    NonStimFR   = PyrPVMeta(i).FiringRate(:,NonStimFlag);
    
    if length(NonStimFR)>4
        PyrPVFR(:,i) = median(NonStimFR,2); %Pyramidal Firing rate is row one, PV is row two
        PyrPVFRMean(:,i) = mean(NonStimFR,2);
    end
end

PyrSSTFR = nan(2,length(PyrSSTMeta));
PyrSSTFRMean = nan(2,length(PyrSSTMeta));
for i = 1:length(PyrSSTMeta)
    NonStimFlag = PyrSSTMeta(i).StimFlag(PyrSSTMeta(i).StimFlag ~= -1) == 0; % The FiringRate of Upstates labeled -1 was not calculated so if the StimFlag has -1 the dimensions won't match with the FR
    NonStimFR   = PyrSSTMeta(i).FiringRate(:,NonStimFlag);
    
    if length(NonStimFR)>4
        PyrSSTFR(:,i) = median(NonStimFR,2); %Pyramidal Firing rate is row one, SST is row two
        PyrSSTFRMean(:,i) = mean(NonStimFR,2);
    end
end

OneTPyrPVFR = nan(2,length(OneTPyrPVMETA));
OneTPyrPVFRMean = nan(2,length(OneTPyrPVMETA));
for i = 1:length(OneTPyrPVMETA)
    NonStimFlag = OneTPyrPVMETA(i).StimFlag(OneTPyrPVMETA(i).StimFlag ~= -1) == 0; % The FiringRate of Upstates labeled -1 was not calculated so if the StimFlag has -1 the dimensions won't match with the FR
    NonStimFR   = OneTPyrPVMETA(i).FiringRate(:,NonStimFlag);
    
    if length(NonStimFR)>4
        OneTPyrPVFR(:,i) = median(NonStimFR,2); %Pyramidal Firing rate is row one, SST is row two
        OneTPyrPVFRMean(:,i) = mean(NonStimFR,2);
    end
end

OneTPyrSSTFR = nan(2,length(OneTPyrSSTMETA));
OneTPyrSSTFRMean = nan(2,length(OneTPyrSSTMETA));
for i = 1:length(OneTPyrSSTMETA)
    NonStimFlag = OneTPyrSSTMETA(i).StimFlag(OneTPyrSSTMETA(i).StimFlag ~= -1) == 0; % The FiringRate of Upstates labeled -1 was not calculated so if the StimFlag has -1 the dimensions won't match with the FR
    NonStimFR   = OneTPyrSSTMETA(i).FiringRate(:,NonStimFlag);
    
    if length(NonStimFR)>4
        OneTPyrSSTFR(:,i) = median(NonStimFR,2); %Pyramidal Firing rate is row one, SST is row two
        OneTPyrSSTFRMean(:,i) = mean(NonStimFR,2);
    end
end

PyrFR    = [PyrPVFR(1,~isnan(PyrPVFR(1,:))) PyrSSTFR(1,~isnan(PyrSSTFR(1,:))) OneTPyrPVFR(1,~isnan(OneTPyrPVFR(1,:))) OneTPyrSSTFR(1,~isnan(OneTPyrSSTFR(1,:)))];
PVFR     = [PyrPVFR(2,~isnan(PyrPVFR(2,:))) OneTPyrPVFR(2,~isnan(OneTPyrPVFR(2,:)))];
SSTFR    = [PyrSSTFR(2,~isnan(PyrSSTFR(2,:))) OneTPyrSSTFR(2,~isnan(OneTPyrSSTFR(2,:)))];

DataFR  = [PyrFR PVFR SSTFR];
IndexFR = [ones(1,length(PyrFR)) 2*ones(1,length(PVFR)) 3*ones(1,length(SSTFR))];

subplot(SP5)
boxplot(DataFR,IndexFR,'Labels',{'PYR','PV','SST'},'Symbol','o')
A = gca;
set(A,'LineWidth',4,'FontSize',16,'FontWeight','bold')
set(A.Children.Children,'Color','k','LineWidth',4)
set(A.Children.Children(1:3),'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1,'MarkerSize',3)
set(A.Children.Children(15:21),'LineStyle','-')
color = {GoodCyan,'r',GoodGreen};
H = findobj(gca,'Tag','Box');
for j=1:length(H)
   patch(get(H(j),'XData'),get(H(j),'YData'),color{j},'FaceAlpha',.5);
end
ylabel('Pyr Up State Duration (s)','FontWeight','bold','FontSize',16)
ylim([-5 40])
box off

[p1,h1]     = ranksum(PyrFR,PVFR);
[p11,h11]   = ranksum(PyrFR,SSTFR);
[p111,h111] = ranksum(PVFR,SSTFR);

if h1==1
    hold on
    line([1.05 1.95],[35 35],'LineWidth',3,'Color','k')
    text(1.4,37,['\ast' '\ast'],'FontWeight','bold','FontSize',20)
end
if h111==1
    hold on
    line([2.05 2.95],[35 35],'LineWidth',3,'Color','k')
    text(2.4,37,['\ast' '\ast'],'FontWeight','bold','FontSize',20)
end
if h11==1
    hold on
    line([1.05 2.95],[38 38],'LineWidth',3,'Color','k')
    text(1.9,40,['\ast' '\ast'],'FontWeight','bold','FontSize',20)
end
%% Panel D
%Getting Correlation Data
PyrPVCorr = nan(1,length(PyrPVMeta));
for i = 1:length(PyrPVMeta)
    NonStimFlag = PyrPVMeta(i).StimFlag(PyrPVMeta(i).StimFlag ~= -1) == 0; % The FiringRate of Upstates labeled -1 was not calculated so if the StimFlag has -1 the dimensions won't match with the FR
    NonStimCorr = PyrPVMeta(i).Corr(:,NonStimFlag);
    
    if length(NonStimCorr)>4
        PyrPVCorr(i) = median(NonStimCorr); %Pyramidal Firing rate is row one, PV is row two
    end
end

OneTPyrPVCorr = nan(1,length(OneTPyrPVMETA));
for i = 1:length(OneTPyrPVMETA)
    NonStimFlag = OneTPyrPVMETA(i).StimFlag(OneTPyrPVMETA(i).StimFlag ~= -1) == 0; % The FiringRate of Upstates labeled -1 was not calculated so if the StimFlag has -1 the dimensions won't match with the FR
    NonStimCorr = OneTPyrPVMETA(i).Corr(:,NonStimFlag);
    
    if length(NonStimCorr)>4
        OneTPyrPVCorr(i) = median(NonStimCorr); %Pyramidal Firing rate is row one, PV is row two
    end
end

PyrSSTCorr = nan(1,length(PyrSSTMeta));
acceptSST = [];
for i = 1:length(PyrSSTMeta)
    NonStimFlag = PyrSSTMeta(i).StimFlag(PyrSSTMeta(i).StimFlag ~= -1) == 0; % The FiringRate of Upstates labeled -1 was not calculated so if the StimFlag has -1 the dimensions won't match with the FR
    NonStimCorr = PyrSSTMeta(i).Corr(:,NonStimFlag);
    
    if length(NonStimCorr)>4
        PyrSSTCorr(i) = median(NonStimCorr); %Pyramidal Firing rate is row one, SST is row two
        acceptSST = [acceptSST i];
    end
end

OneTPyrSSTCorr = nan(1,length(OneTPyrSSTMETA));
for i = 1:length(OneTPyrSSTMETA)
    NonStimFlag = OneTPyrSSTMETA(i).StimFlag(OneTPyrSSTMETA(i).StimFlag ~= -1) == 0; % The FiringRate of Upstates labeled -1 was not calculated so if the StimFlag has -1 the dimensions won't match with the FR
    NonStimCorr = OneTPyrSSTMETA(i).Corr(:,NonStimFlag);
    
    if length(NonStimCorr)>4
        OneTPyrSSTCorr(i) = median(NonStimCorr); %Pyramidal Firing rate is row one, PV is row two
    end
end

PyrPVCORR     = [PyrPVCorr(~isnan(PyrPVCorr)) OneTPyrPVCorr(~isnan(OneTPyrPVCorr))];
PyrSSTCORR    = [PyrSSTCorr(~isnan(PyrSSTCorr)) OneTPyrSSTCorr(~isnan(OneTPyrSSTCorr))];

DataCorr  = [PyrPVCORR PyrSSTCORR];
IndexCorr = [1*ones(1,length(PyrPVCORR)) 2*ones(1,length(PyrSSTCORR))];

[p2,h2,stats2]     = ranksum(PyrPVCORR,PyrSSTCORR);

subplot(SP6)
boxplot(DataCorr,IndexCorr,'Labels',{'PYR-PV','PYR-SST'},'Symbol','o')
A = gca;
set(A,'LineWidth',4,'FontSize',16,'FontWeight','bold')
set(A.Children.Children,'Color','k','LineWidth',4)
set(A.Children.Children(1:2),'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1,'MarkerSize',3)
set(A.Children.Children(11:14),'LineStyle','-')
color = {GoodCyan,'r'};
H = findobj(gca,'Tag','Box');
for j=1:length(H)
   patch(get(H(j),'XData'),get(H(j),'YData'),color{j},'FaceAlpha',.5);
end
ylabel('Correlation','FontWeight','bold','FontSize',16)
ylim([0 1])
box off


fprintf('Mean of Median FR Pyr/PV = %4.2f/%4.2f\n',nanmean([PyrPVFR OneTPyrPVFR] ,2));
fprintf('Mean of Median FR Pyr/SST = %4.2f/%4.2f\n',nanmean([PyrSSTFR OneTPyrSSTFR],2));
fprintf('Mean of Merged Median FR Pyr = %4.2f\n',nanmean(PyrFR,2));
fprintf('N of FR Pyr = %4.2f\n',length(PyrFR));
fprintf('N of FR PV = %4.2f\n',length(PVFR));
fprintf('N of FR SST %4.2f\n',length(SSTFR));








