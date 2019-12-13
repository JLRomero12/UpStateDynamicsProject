clear
close all

cd('C:\Users\BuonoLab\Desktop\Juan Luis\PV_SST_Opto_Study_19\Misc-UnOrganizedMess\SimultaneousRecording')
load META

ExampleCell = 28;
ExampleCond = 1;
ExampleTrace = 6;

%% INIT GRAPHICS
MainFigure = figure('Name','MainFigure');
set(gcf,'color','w');
ScreenSize = get(groot,'ScreenSize');
MainFigure.Position = [0 0 ScreenSize(3)/1.2 ScreenSize(4)/1.2]; %This is to Force the figure to have a vertical orientation

SP1 = subplot('Position',[0.02,0.52,0.73,0.46]);
SP11 = subplot('Position',[0.02,0.01,0.73,0.46]);
SP2 = subplot('Position',[0.81,0.55,0.18,0.4]);
SP3 = subplot('Position',[0.81,0.03,0.18,0.4]);

% INSET2 = axes('Position',[0.85,0.7,0.145,0.28]);
INSET1 = axes('Position',[0.01,0.75,0.14,0.21]);
INSET2 = axes('Position',[0.01,0.26,0.14,0.21]);

%% Example Upstates
Files = dir('C*.mat');
numFiles = length(Files);
load(Files(ExampleCell).name)

dt = SampleRate(ExampleCond);
TracePY = CH1{ExampleCond}(ExampleTrace,60000:67000);
TracePV = CH2{ExampleCond}(ExampleTrace,60000:67000);
X = (1:length(TracePY))*(dt/1000);

subplot(SP1)
plot(X,TracePY,'k','LineWidth',1.25)
hold on
plot(X,TracePV,'r','LineWidth',1.25)
rectangle('Position',[0.22 -65 0.18 29],'FaceColor','none','EdgeColor','k','LineStyle','--','LineWidth',2.0)
xlim([0 1.41])
ylim([-67.5 5])
text(0,-67,'-65mV','FontWeight','bold')
line([1.21 1.41],[-65 -65],'Color','k','LineWidth',1.25)
text(1.3,-66    ,'200ms','FontWeight','bold')
line([1.41 1.41],[-65 -45],'Color','k','LineWidth',1.25)
text(1.424,-51,'20mV','FontWeight','bold')
A = get(gca,'Children');
A(1).Rotation = 270;
hold off
axis off
box off

subplot(SP11)
hold on
rectangle('Position',[0.22 -65 0.18 29],'FaceColor','none','EdgeColor','k','LineStyle','--','LineWidth',2.0)
xlim([0 1.41])
ylim([-66.5 15])
text(0,-66,'-65mV','FontWeight','bold')
line([1.21 1.41],[-64 -64],'Color','k','LineWidth',1.25)
text(1.3,-65,'200ms','FontWeight','bold')
line([1.41 1.41],[-64 -44],'Color','k','LineWidth',1.25)
text(1.424,-51,'20mV','FontWeight','bold')
A = get(gca,'Children');
A(1).Rotation = 270;
hold off
axis off
box off

axes(INSET1)

plot(X(1000:2000),TracePY(1000:2000),'k','LineWidth',1.25)
hold on
plot(X(1000:2000),TracePV(1000:2000),'r','LineWidth',1.25)
scatter([X(1370) X(1312)],[TracePY(1370) TracePV(1312)],20,[0 0 0;1 0 0],'filled')
line([X(1370) X(1370)],[-67 TracePY(1370)],'Color',[0 0 1 0.5],'LineStyle','-','LineWidth',1.25)
line([X(1000) X(1370)],[TracePY(1370) TracePY(1370)],'Color',[0 0 1 0.5],'LineStyle','-','LineWidth',1.25)
line([X(1312) X(1312)],[-67 TracePV(1312)],'Color',[0 0 1 0.5],'LineStyle','-','LineWidth',1.25)
line([X(1000) X(1312)],[TracePV(1312) TracePV(1312)],'Color',[0 0 1 0.5],'LineStyle','-','LineWidth',1.25)
hold off
xlim([0.22 0.4])
ylim([-66 -36])
MainFigure.Children(1).YTick = [];
MainFigure.Children(1).XTick = [];

%% Average UpStart Diff
eventind = vertcat(META.numEvents)>9; %More than 9 events

AllStartDiff = [META(eventind).MedianStartDiff];
[h,p,s,stat] = ttest(AllStartDiff);

% subplot(SP2)
axes(SP2)
boxplot(AllStartDiff)
ylabel('Onset PYR - Onset PV','FontWeight','bold','FontSize',12)
A = gca;
set(A,'LineWidth',2.0,'FontSize',10,'FontWeight','bold','XTick',[])
set(A.Children.Children,'LineWidth',2.0,'Color','k')
A = get(gca,'Children');
A.Children(2).XData = [0.8 1.2];
A.Children(3).XData = [0.8 0.8 1.2 1.2 0.8];
color = 'k';
H = findobj(gca,'Tag','Box');
for j=1:length(H)
    patch(get(H(j),'XData'),get(H(j),'YData'),color(j),'FaceAlpha',.5);
end
% A = get(gca,'Children');
% A.Children(2).XData = [0.7 1.3];
% A.Children(3).XData = [0.7 0.7 1.3 1.3 0.7];
title('Difference UpStart PV','FontSize',14)
ylim([-90 250])
box off

if h==1
    hold on
    text(0.97,223,'\ast','FontWeight','bold','FontSize',12)
end

% %% XCorr PeakShift
% PS = [META(eventind).MeanPeakShifts];
% PS50 = PS(1,:);
% PS100 = PS(2,:);
% 
% [h50,p50] = ttest(PS50);
% [h100,p100,s100,stat100] = ttest(PS100);
% 
% subplot(SP3)
% boxplot([PS100'],'Labels',{'XCorr'})
% ylabel('XCorrelation Peak','FontWeight','bold','FontSize',12)
% A = gca;
% set(A,'LineWidth',2.0,'FontSize',12,'FontWeight','bold')
% set(A.Children.Children,'LineWidth',2.0,'Color','k')
% A = get(gca,'Children');
% A.Children(2).XData = [0.8 1.2];
% A.Children(3).XData = [0.8 0.8 1.2 1.2 0.8];
% color = 'k';
% H = findobj(gca,'Tag','Box');
% for j=1:length(H)
%     patch(get(H(j),'XData'),get(H(j),'YData'),color(j),'FaceAlpha',.5);
% end
% title('XCorrelation 100ms','FontSize',14)
% ylim([-5 15])
% box off
% 
% if h100==1
%     hold on
%     text(0.85,13.5,'p<0.05','FontWeight','bold','FontSize',12)
% end
% 
% print -dtiffn Figure4A

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% SOM DIff in UpState %%%%%%%%%%%%%%%%%%%%%%%%%
% % INIT GRAPHICS
% SOMFigure = figure('Name','MainFigure');
% set(gcf,'color','w');
% ScreenSize = get(groot,'ScreenSize');
% SOMFigure.Position = [0 0 ScreenSize(3)/1.2 ScreenSize(4)/1.2]; %This is to Force the figure to have a vertical orientation
% 
% SP1 = subplot('Position',[0.02,0.5,0.73,0.48]);
% SP11 = subplot('Position',[0.02,0.03,0.73,0.46]);
% SP2 = subplot('Position',[0.81,0.55,0.18,0.4]);
% SP3 = subplot('Position',[0.81,0.08,0.18,0.4]);
% 
% % INSET2 = axes('Position',[0.85,0.7,0.145,0.28]);
% INSET1 = axes('Position',[0.01,0.75,0.14,0.21]);

%% Getting Data
cd('C:\Users\BuonoLab\Desktop\Juan Luis\MainStudy\SOM-tdTomato')

load SOMTDTomatoMasterMETA
for i = 1:length(MasterMETA)
    if i~=5 %Cell 5 not Dual Recording
        tracePYR = squeeze(MasterMETA(i).UpStates(1,:,:));
        traceSOM = squeeze(MasterMETA(i).UpStates(2,:,:));
        STARTDIFF = [];
        
        for j = 1:size(tracePYR,1)
            thrPYR = mode(tracePYR(j,:))+Param.thr;
            thrSOM = mode(traceSOM(j,:))+Param.thr;
            window = 1:Param.preEventTime*2;
            MaxDiff = Param.MaxIntEvent; %30 ms
            
            AbovePY = find(tracePYR(j,window)>(thrPYR));
            AboveSOM = find(traceSOM(j,window)>(thrSOM));

            if ~isempty(AbovePY) && ~isempty(AboveSOM)
                OnsetStartPY = (find(diff([-9999 AbovePY])>MaxDiff));
                OffsetStartPY=AbovePY([OnsetStartPY(2:length(OnsetStartPY))-1 length(AbovePY)]);
                OnsetStartPY = AbovePY(OnsetStartPY);
                [~,n] = max(OffsetStartPY-OnsetStartPY);
                UpStartPY = OnsetStartPY(n);
                clear n
                
                OnsetStartSOM = (find(diff([-9999 AboveSOM])>MaxDiff));
                OffsetStartSOM=AboveSOM([OnsetStartSOM(2:length(OnsetStartSOM))-1 length(AboveSOM)]);
                OnsetStartSOM = AboveSOM(OnsetStartSOM);
                [~,n] = max(OffsetStartSOM-OnsetStartSOM);
                UpStartSOM = OnsetStartSOM(n);
                
                startdiff = UpStartPY-UpStartSOM;
                STARTDIFF = [STARTDIFF startdiff];
                
            end
    
        end
    end
    MasterMETA(i).StartDiff = STARTDIFF;
    MasterMETA(i).MedianStartDiff = median(STARTDIFF);

end

tracePYR = squeeze(MasterMETA(2).UpStates(1,:,:));
traceSOM = squeeze(MasterMETA(2).UpStates(2,:,:));

ExamplePYR = [tracePYR(4,5000:5500) tracePYR(2,1:3680)];
ExampleSOM = [(traceSOM(4,5000:5500)-1) traceSOM(2,1:3680)];

X = (1:length(ExamplePYR))*(Param.dt/1000);

subplot(SP11)
plot(X,ExamplePYR,'k','LineWidth',1.25)
hold on
plot(X,ExampleSOM,'b','LineWidth',1.25)
rectangle('Position',[0.26 -65 0.18 29],'FaceColor','none','EdgeColor','k','LineStyle','--','LineWidth',2.0)
xlim([0 1.75])
ylim([-67.5 5])
text(0,-64,'-62mV','FontWeight','bold')
hold off
axis off
box off

line([1.56 1.75],[-65 -65],'Color','k','LineWidth',1.25)
text(1.615,-66,'200ms','FontWeight','bold')
line([1.75 1.75],[-65 -45],'Color','k','LineWidth',1.25)
text(1.765,-51,'20mV','FontWeight','bold')
A = get(gca,'Children');
A(1).Rotation = 270;

axes(INSET2)
plot(X(520:880),ExamplePYR(520:880),'k','LineWidth',1.25)
hold on
plot(X(520:880),ExampleSOM(520:880),'b','LineWidth',1.25)
scatter([X(602) X(605)],[ExamplePYR(602) ExampleSOM(605)],20,[0 0 0;0 0 1],'filled')
line([X(602) X(602)],[-64 ExamplePYR(602)],'Color',[1 0 0 0.5],'LineStyle','-','LineWidth',1.25)
line([X(520) X(602)],[ExamplePYR(602) ExamplePYR(602)],'Color',[1 0 0 0.5],'LineStyle','-','LineWidth',1.25)
line([X(605) X(605)],[-64 ExampleSOM(605)],'Color',[1 0 0 0.5],'LineStyle','-','LineWidth',1.25)
line([X(520) X(605)],[ExampleSOM(605) ExampleSOM(605)],'Color',[1 0 0 0.5],'LineStyle','-','LineWidth',1.25)
hold off
xlim([0.26 0.44])
ylim([-64 -32])
MainFigure.Children(1).YTick = [];
MainFigure.Children(1).XTick = [];


%% Difference in UpStart
AllStartDiff = [MasterMETA(:).MedianStartDiff];
[h,p,s,stat] = ttest(AllStartDiff);

% subplot(SP2)
axes(SP3)
boxplot(AllStartDiff)
ylabel('Onset PYR - Onset SST','FontWeight','bold','FontSize',12)
A = gca;
set(A,'LineWidth',2.0,'FontSize',10,'FontWeight','bold','XTick',[])
set(A.Children.Children,'LineWidth',2.0,'Color','k')
A = get(gca,'Children');
A.Children(2).XData = [0.8 1.2];
A.Children(3).XData = [0.8 0.8 1.2 1.2 0.8];
color = 'k';
H = findobj(gca,'Tag','Box');
for j=1:length(H)
    patch(get(H(j),'XData'),get(H(j),'YData'),color(j),'FaceAlpha',.5);
end
title('Difference UpStart SST','FontSize',14)
ylim([-65 5])
box off

if h==1
    hold on
    text(0.85,223,'\ast','FontWeight','bold','FontSize',12)
end

print -dtiffn Figure4
