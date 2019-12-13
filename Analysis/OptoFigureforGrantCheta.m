close all
clear

cd('C:\Users\BuonoLab\Desktop\Juan Luis\MainStudy\OptogeneticFigureforGRant')

load('PVChetaMETA.mat')

%Variables for the Duration of UpState
PVExampleCell = 4;
PVStimUpStates = [17,21,28];
PVNonUpStates = [34,37,50];
SOMExample = 0;
MinDuration = 0.25;

%Variables for the Difference in Voltage
dt = Param.dt;
chunksize = 275/dt;
minlength = 250/dt;
stimdelay = 250/dt;

PVMETA = MasterMETA;

%% INIT GRAPHICS
MainFigure = figure('Name','MainFigure');
set(gcf,'color','w');
ScreenSize = get(groot,'ScreenSize');
MainFigure.Position = [50 100 ScreenSize(3)/1.8 ScreenSize(4)/1.3]; %This is to Force the figure to have a vertical orientation

SP1 = subplot('Position',[0.01 0.7 0.95 0.3]);
SP2 = subplot('Position',[0.01 0.4 0.95 0.3]);
SP3 = subplot('Position',[0.1 0.05 0.3 0.3]);
SP4 = subplot('Position',[0.6 0.05 0.3 0.3]);

colormap(jet)
%% Panel A
StimInd = find(MasterMETA(PVExampleCell).StimFlag==1); StimInd = StimInd(PVStimUpStates);
StimUp = squeeze(MasterMETA(PVExampleCell).UpStates(1,StimInd,:));
PreUpState(1,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,StimInd(1)-1,12750:13650))';
PreUpState(2,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,StimInd(1)-1,7270:8170))';
PreUpState(3,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,StimInd(1)-1,6530:7430))';
ISI = 0;

subplot(SP2)
hold on
plot([PreUpState(1,:) StimUp(1,1:5500)],'Color',[0.1 0.1 0.1],'LineWidth',2.0)
plot([PreUpState(2,:) StimUp(2,1:5500)],'Color',[0.35 0.35 0.35],'LineWidth',2.0)
plot([PreUpState(3,:) StimUp(3,1:5500)],'Color',[0.6 0.6 0.6],'LineWidth',2.0)

%rectangle('Position',[600 -72 1000 60],'FaceColor',[0.5 0.7 1 0.5],'EdgeColor','none');
for k = 1:25
    rectangle('Position',[1500+ISI -72 20 29],'FaceColor',[0.5 0.7 1 0.5],'EdgeColor','none');
    ISI = ISI+40;
end
hold off
xlim([0 6400])
ylim([-72 -35])
box off
axis off


%% Panel B
NonStimInd = find(MasterMETA(PVExampleCell).StimFlag==0); NonStimInd = NonStimInd(PVNonUpStates);
NonStimUp = squeeze(MasterMETA(PVExampleCell).UpStates(1,NonStimInd,:));
PreUpState(1,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,NonStimInd(1)-1,2800:3700))';
PreUpState(2,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,NonStimInd(2)-1,950:1850))';
PreUpState(3,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,NonStimInd(2)-1,6450:7350))';


subplot(SP1)
hold on
plot([PreUpState(1,:) NonStimUp(1,1:5500)],'Color',[0.1 0.1 0.1],'LineWidth',2.0)
plot([PreUpState(2,:) NonStimUp(2,1:5500)],'Color',[0.35 0.35 0.35],'LineWidth',2.0)
plot([PreUpState(3,:) NonStimUp(3,1:5500)],'Color',[0.6 0.6 0.6],'LineWidth',2.0)
rectangle('Position',[1500 -72 1000 29],'FaceColor','none','EdgeColor','k','LineStyle','--','LineWidth',2.0);
line([5300 6300],[-48 -48],'Color','k','LineWidth',2.0)
line([6300 6300],[-48 -38],'Color','k','LineWidth',2.0)
text(5575,-49,'500ms','FontWeight','bold','FontSize',12)
text(6400,-40,'10mV','FontWeight','bold','FontSize',12,'Rotation',270)
hold off
xlim([0 6400])
ylim([-72 -35])
box off
axis off

%% Panel C
for i = 1:length(MasterMETA) %Organazing the Durations of the UpStates
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
    
    if length(MasterMETA(i).UpDurPY)>9 && i~=10 && i~=11 %Finding the Volage during stimmed and non stimmed
        CH1 = squeeze(MasterMETA(i).UpStates(1,:,:));
        CH3 = squeeze(MasterMETA(i).UpStates(3,:,:));
        UpDurPY = MasterMETA(i).UpDurPY;
        chunk = zeros(size(CH3,1),chunksize*3.5);
        
        %% Getting chunk centered around the first stimulus or 250 ms from onset (in non stimulated)
        for j = 1:size(CH3,1)
            if UpDurPY(j)>minlength
                ind = find(CH3(j,1:chunksize*3)>10,1);
                if ~isempty(ind) && ind>500
                    chunk(j,:) = CH1(j,ind-chunksize:ind+(chunksize*2.5)-1);
                elseif isempty(ind)
                    chunk(j,:) = CH1(j,1+((Param.preEventTime+stimdelay)-chunksize):(Param.preEventTime+stimdelay)+(chunksize*2.5));
                end
            end
        end
        
        stimchunk        = chunk(StimIndex,:);
        stimchunk        = stimchunk(any(stimchunk,2),:);
        filtstimchunk    = medfilt1(stimchunk(:,:)',5/dt,'truncate')';
        StimOnVolt       = mean(filtstimchunk(:,chunksize-(10/dt):chunksize),2);
        StimOffVolt      = mean(filtstimchunk(:,chunksize+(500/dt):chunksize+(510/dt)),2);
        
        nonstimchunk     = chunk(NonStimIndex,:);
        nonstimchunk     = nonstimchunk(any(nonstimchunk,2),:);
        filtnonstimchunk = medfilt1(nonstimchunk(:,:)',5/dt,'truncate')';
        NonStimOnVolt    = mean(filtnonstimchunk(:,chunksize-(10/dt):chunksize),2);
        NonStimOffVolt   = mean(filtnonstimchunk(:,chunksize+(500/dt):chunksize+(510/dt)),2);
        
        PVMETA(i).OnVolt  = [mean(StimOnVolt); mean(NonStimOnVolt)];
        PVMETA(i).OffVolt = [mean(StimOffVolt); mean(NonStimOffVolt)];
    end
end

NonDur  = [MasterMETA(:).PassNonStimUpDur];
StimDur = [MasterMETA(:).PassStimUpDur];

[p1,h1] = signrank(NonDur,StimDur);

subplot(SP3)
boxplot([NonDur' StimDur'],'Labels',{'Off','On'},'Symbol','o')
title('PV','FontSize',8)
ylim([0 6.2])
A = gca;
set(A,'LineWidth',2,'FontSize',12,'FontWeight','bold')
set(A.Children.Children,'Color','k','LineWidth',2)
set(A.Children.Children(1:2),'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1,'MarkerSize',3)
color = {[0.5 0.7 1],'k'};
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
SOMMETA = MasterMETA;

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
    
    if length(MasterMETA(i).UpDurPY)>9 && i~=10 && i~=11 %Finding the Volage during stimmed and non stimmed
        CH1 = squeeze(MasterMETA(i).UpStates(1,:,:));
        CH3 = squeeze(MasterMETA(i).UpStates(3,:,:));
        UpDurPY = MasterMETA(i).UpDurPY;
        chunk = zeros(size(CH3,1),chunksize*3.5);
        
        %% Getting chunk centered around the first stimulus or 250 ms from onset (in non stimulated)
        for j = 1:size(CH3,1)
            if UpDurPY(j)>minlength
                ind = find(CH3(j,1:chunksize*3)>10,1);
                if ~isempty(ind) && ind>500
                    chunk(j,:) = CH1(j,ind-chunksize:ind+(chunksize*2.5)-1);
                elseif isempty(ind)
                    chunk(j,:) = CH1(j,1+((Param.preEventTime+stimdelay)-chunksize):(Param.preEventTime+stimdelay)+(chunksize*2.5));
                end
            end
        end
        
        stimchunk        = chunk(StimIndex,:);
        stimchunk        = stimchunk(any(stimchunk,2),:);
        filtstimchunk    = medfilt1(stimchunk(:,:)',5/dt,'truncate')';
        StimOnVolt       = mean(filtstimchunk(:,chunksize-(10/dt):chunksize),2);
        StimOffVolt      = mean(filtstimchunk(:,chunksize+(500/dt):chunksize+(510/dt)),2);
        
        nonstimchunk     = chunk(NonStimIndex,:);
        nonstimchunk     = nonstimchunk(any(nonstimchunk,2),:);
        filtnonstimchunk = medfilt1(nonstimchunk(:,:)',5/dt,'truncate')';
        NonStimOnVolt    = mean(filtnonstimchunk(:,chunksize-(10/dt):chunksize),2);
        NonStimOffVolt   = mean(filtnonstimchunk(:,chunksize+(500/dt):chunksize+(510/dt)),2);
        
        SOMMETA(i).OnVolt  = [mean(StimOnVolt); mean(NonStimOnVolt)];
        SOMMETA(i).OffVolt = [mean(StimOffVolt); mean(NonStimOffVolt)];
    end
end

NonDur = [MasterMETA(:).PassNonStimUpDur];
StimDur = [MasterMETA(:).PassStimUpDur];

[p2,h2] = signrank(NonDur,StimDur);

subplot(SP4)
boxplot([NonDur' StimDur'],'Labels',{'Off','On'},'Symbol','o')
title('SOM','FontSize',8)
ylim([0 6.2])
A = gca;
set(A,'LineWidth',2,'FontSize',12,'FontWeight','bold')
set(A.Children.Children,'Color','k','LineWidth',2)
set(A.Children.Children(1:2),'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1,'MarkerSize',3)
color = {[0.5 0.7 1],'k'};
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

print '-dtiffn' OptoFigureAB

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Change in Voltage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphics
PanelC = figure('Name','MainFigure');
set(gcf,'color','w');
ScreenSize = get(groot,'ScreenSize');
PanelC.Position = [50 100 ScreenSize(3)/1.8 ScreenSize(4)/1.3]; %This is to Force the figure to have a vertical orientation

SP1 = subplot('Position',[0.05 0.05 0.4 0.3]);
SP2 = subplot('Position',[0.55 0.05 0.4 0.3]);


%% Getting Data
OnPVVoltage   = [PVMETA(:).OnVolt];
OffPVVoltage  = [PVMETA(:).OffVolt];
OnSOMVoltage  = [SOMMETA(:).OnVolt];
OffSOMVoltage = [SOMMETA(:).OffVolt];

subplot(SP1)
boxplot([OnPVVoltage(2,:)' OnPVVoltage(1,:)' OffPVVoltage(2,:)' OffPVVoltage(1,:)'],'Labels',{'PreNon','PreStim','PostNon','PostStim'},'Symbol','o')
title('Difference in Voltage PV ChETA','FontSize',12)
ylabel('Voltage (mV)')
A = gca;
set(A,'LineWidth',2,'FontSize',12,'FontWeight','bold')
set(A.Children.Children,'Color','k','LineWidth',2)
set(A.Children.Children(3:4),'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1,'MarkerSize',3)
color = {[0.5 0.7 1],'k',[0.5 0.7 1],'k'};
H = findobj(gca,'Tag','Box');
for j=1:length(H)
   patch(get(H(j),'XData'),get(H(j),'YData'),color{j},'FaceAlpha',.5);
end
ylim([-70 -43])
box off

[h1,p1] = ttest(OnPVVoltage(2,:),OnPVVoltage(1,:));
[h11,p11] = ttest(OffPVVoltage(2,:),OffPVVoltage(1,:));

if h11 == 1
    hold on
    line([3 4],[-45 -45],'Color','k','LineWidth',2) 
    text(3.32,-44,'\ast \ast','FontWeight','bold','FontSize',14)
end

subplot(SP2)
boxplot([OnSOMVoltage(2,:)' OnSOMVoltage(1,:)' OffSOMVoltage(2,:)' OffSOMVoltage(1,:)'],'Labels',{'PreNon','PreStim','PostNon','PostStim'})
title('Difference in Voltage SOM-ChETA','FontSize',12)
ylabel('Voltage (mV)')
A = gca;
set(A,'LineWidth',2,'FontSize',12,'FontWeight','bold')
set(A.Children.Children,'Color','k','LineWidth',2)
set(A.Children.Children(1:2),'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1)
color = {[0.5 0.7 1],'k',[0.5 0.7 1],'k'};
H = findobj(gca,'Tag','Box');
for j=1:length(H)
   patch(get(H(j),'XData'),get(H(j),'YData'),color{j},'FaceAlpha',.5);
end
ylim([-70 -43])
box off

[h2,p2] = ttest(OnSOMVoltage(2,:),OnSOMVoltage(1,:));
[h22,p22] = ttest(OffSOMVoltage(2,:),OffSOMVoltage(1,:));

if h22 == 1
    hold on
    line([3 4],[-45 -45],'Color','k','LineWidth',2) 
    text(3.2,-43,'\ast \ast','FontWeight','bold','FontSize',14)
end

print '-dtiffn' OptoFigureC
