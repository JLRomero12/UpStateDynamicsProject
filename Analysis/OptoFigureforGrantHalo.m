close all
clear

cd('C:\Users\BuonoLab\Desktop\Juan Luis\MainStudy\OptogeneticFigureforGRant')
load('PVHaloMETA.mat')

MinDuration = 0.25;

%% INIT GRAPHICS
MainFigure = figure('Name','MainFigure');
set(gcf,'color','w');
ScreenSize = get(groot,'ScreenSize');
MainFigure.Position = [50 100 ScreenSize(3)/1.7 ScreenSize(4)/2]; %This is to Force the figure to have a vertical orientation

SP1 = subplot('Position',[0.1 0.05 0.3 0.8]);
SP2 = subplot('Position',[0.6 0.05 0.3 0.8]);

%% Panel C
for i = 1:length(MasterMETA)
    if length(MasterMETA(i).UpDurPY)>9
        StimIndex       = MasterMETA(i).StimFlag==1;
        NonStimIndex    = MasterMETA(i).StimFlag==0;
        
        NonStimUpDur    = (MasterMETA(i).UpDurPY(1,NonStimIndex,:))*(Param.dt/1000);
        StimUpDur       = (MasterMETA(i).UpDurPY(1,StimIndex,:))*(Param.dt/1000);
        PassNon         = NonStimUpDur(NonStimUpDur>MinDuration);
        PassStim        = StimUpDur(StimUpDur>MinDuration);
    end
    MasterMETA(i).PassNonStimUpDur = median(PassNon);
    MasterMETA(i).PassStimUpDur = median(PassStim);
end

NonDur = [MasterMETA(:).PassNonStimUpDur];
StimDur = [MasterMETA(:).PassStimUpDur];

[p1,h1] = signrank(NonDur,StimDur);

subplot(SP1)
boxplot([NonDur' StimDur'],'Labels',{'Off','On'},'Symbol','o')
title('PV','FontSize',8)
ylim([0 6.2])
A = gca;
set(A,'LineWidth',2,'FontSize',12,'FontWeight','bold')
set(A.Children.Children,'Color','k','LineWidth',2)
set(A.Children.Children(1:2),'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1,'MarkerSize',3)
color = {[1 0.5 0.1],'k'};
H = findobj(gca,'Tag','Box');
for j=1:length(H)
   patch(get(H(j),'XData'),get(H(j),'YData'),color{j},'FaceAlpha',.5);
end
ylabel('Pyr Up State Duration (s)','FontWeight','bold','FontSize',12)
box off

if h1==1
    hold on
    line([1 2],[5 5],'LineWidth',2,'Color','k')
    text(1.4,5.3,['\ast' '\ast'],'FontWeight','bold','FontSize',16)
end

%% Panel D

clear MasterMETA
load('SOMChetaMETA.mat')

for i = 1:length(MasterMETA)
    if length(MasterMETA(i).UpDurPY)>9
        StimIndex       = MasterMETA(i).StimFlag==1;
        NonStimIndex    = MasterMETA(i).StimFlag==0;
        
        NonStimUpDur    = (MasterMETA(i).UpDurPY(1,NonStimIndex,:))*(Param.dt/1000);
        StimUpDur       = (MasterMETA(i).UpDurPY(1,StimIndex,:))*(Param.dt/1000);
        PassNon         = NonStimUpDur(NonStimUpDur>MinDuration);
        PassStim        = StimUpDur(StimUpDur>MinDuration);
    end
    MasterMETA(i).PassNonStimUpDur = median(PassNon);
    MasterMETA(i).PassStimUpDur = median(PassStim);
end

NonDur = [MasterMETA(:).PassNonStimUpDur];
StimDur = [MasterMETA(:).PassStimUpDur];

[p2,h2] = signrank(NonDur,StimDur);

subplot(SP2)
boxplot([NonDur' StimDur'],'Labels',{'Off','On'},'Symbol','o')
title('SOM','FontSize',8)
ylim([0 6.2])
A = gca;
set(A,'LineWidth',2,'FontSize',12,'FontWeight','bold')
set(A.Children.Children,'Color','k','LineWidth',2)
set(A.Children.Children(1:2),'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1,'MarkerSize',3)
color = {[1 0.5 0.1],'k'};
H = findobj(gca,'Tag','Box');
for j=1:length(H)
   patch(get(H(j),'XData'),get(H(j),'YData'),color{j},'FaceAlpha',.5);
end
ylabel('Pyr Up State Duration (s)','FontWeight','bold','FontSize',12)
box off

if h2==1
    hold on
    line([1 2],[5 5],'LineWidth',2,'Color','k')
    text(1.4,5.3,['\ast' '\ast'],'FontWeight','bold','FontSize',16)
end
