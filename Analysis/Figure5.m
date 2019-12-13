% clear
close all

load MasterMETA
META = MasterMETA
dt = 0.5;


%% Bullian Indexes & Extracting Data
%Allocating
transind = zeros(length(META),1);
StimFlag = cell(length(META),1);
FiringRate = cell(length(META),1);
StimUpDur = cell(length(META),1);
NonStimUpDur= cell(length(META),1);

%More than 9 events
eventind = vertcat(META.numEvents)>9; 

%Transfection type is flex-tdTomato DIO-Ch2r
for i = 1:length(META) 
    transind(i,1) = strcmp(META(i).TransType,'DIO-Cheta');
    StimFlag{i} = META(i).StimFlag;
    FiringRate{i} = META(i).FiringRate;
    StimUpDur{i} = META(i).StimUpDurPY;
    NonStimUpDur{i} = META(i).NonStimUpDurPY;
end

StimFlag = StimFlag(eventind&transind);
FiringRate = FiringRate(eventind&transind);
StimUpDur = StimUpDur(eventind&transind);
NonStimUpDur = NonStimUpDur(eventind&transind);
AllFR = zeros(4,sum(eventind&transind));

for j = 1:sum(eventind&transind)
    AllFR(1:2,j) = median(FiringRate{j}(:,~StimFlag{j}),2); %FR of entire unstimmulated upstates row1 PYR and row2 PV
    AllFR(3:4,j) = median(FiringRate{j}(:,StimFlag{j}),2); %FR of entire stimmulated upstates row3 PYR and row4 PV
end
%% PANEL E Duration of Stimmed Upstate vs Non-Stimmed Upstate

StimUpDur = (cellfun(@median,StimUpDur)).*(dt/1000);
NonStimUpDur = (cellfun(@median,NonStimUpDur)).*(dt/1000);

figure
boxplot([NonStimUpDur StimUpDur],'Labels',{'Non-Stim','Stim'},'Symbol','o')
title('Avrg Duration Stim vs Non-Stim','FontSize',12)
A = gca;
set(A,'LineWidth',2,'FontSize',12,'FontWeight','bold')
set(A.Children.Children,'Color','k','LineWidth',2)
set(A.Children.Children(1:2),'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1)
color = ['k','k'];
H = findobj(gca,'Tag','Box');
for j=1:length(H)
   patch(get(H(j),'XData'),get(H(j),'YData'),color(j),'FaceAlpha',.5);
end
ylabel('PYR Up State Duration (s)','FontWeight','bold','FontSize',12)
box off

% figure
% boxplot([StimUpDur NonStimUpDur],'Labels',{'Duration Stim','Duration Non-Stim'})
% [h4,p4,c4,stats4] = ttest2(StimUpDur',NonStimUpDur');
% title(['Average Duration of UpState p=' num2str(p4)])

% figure
% deltaDur = NonStimUpDur-StimUpDur;
% deltaPVFR = AllFR(4,:)-AllFR(3,:);
% x = deltaPVFR;
% y=deltaDur;
% scatter(deltaPVFR,deltaDur)
% hold on
% [B,BINT,R,RINT,STATS] = regress(y,[ones(1,length(x)); x]');
% plot(x,B(1)+x*B(2),'k','linewidth',2)
