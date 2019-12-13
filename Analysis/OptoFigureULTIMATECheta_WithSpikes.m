close all
clear
warning('off','signal:findpeaks:largeMinPeakHeight')

cd('C:\Users\IzquierdoLab\Desktop\Juan Luis\PV_SST_Opto_Study_19\OptoUpStateData\PV-Cheta\DoubleTrans\')
load('PVCheta2TransMETA.mat')

%Variables for the Duration of UpState
dt             = Param.dt;
PVExampleCell  = 2;
PVStimUpStates = [1,2,4];
PVNonUpStates  = [3,5,7];
PVPosCtrl      = 'C:\Users\IzquierdoLab\Desktop\Juan Luis\PV_SST_Opto_Study_19\OptoUpStateData\PV-Cheta\DoubleTrans\PVOnlyS5C3_0730_.mat';
SSTPosCtrl     = 'C:\Users\IzquierdoLab\Desktop\Juan Luis\PV_SST_Opto_Study_19\OptoUpStateData\SOM-Cheta\DoubleTrans\S4C1_0519_SOMOnly.mat';
MinDuration    = 0.25;
FailMinLen     = 0.750;

%Variables for the Difference in Voltage
dt = Param.dt;
stimlen = 500/dt;
prestimsize = 275/dt;
poststimsize = 5100/dt;
startlight = prestimsize;
minlength = 250/dt;
stimdelay = 250/dt;
voltwinsize = 100/dt;
filtorder = 10/dt;
voltstimlen = 215/dt;
SSTGroup = 1:30;
% SSTGroup = 21:30;

%Variables for Firing Rate
spikethr = -25;
minlengthFR = 250/dt;


PVMETA = MasterMETA;

%% INIT GRAPHICS
MainFigure = figure('Name','MainFigure');
set(gcf,'color','w');
ScreenSize = get(groot,'ScreenSize');
MainFigure.Position = [50 30 ScreenSize(3)/2 ScreenSize(4)/1.1]; %This is to Force the figure to have a vertical orientation

SP1 = subplot('Position',[0.01 0.81 0.95 0.21]);
SP2 = subplot('Position',[0.01 0.59 0.95 0.21]);
SP3 = subplot('Position',[0.01 0.33 0.24 0.22]);
SP4 = subplot('Position',[0.51 0.33 0.24 0.22]);
SP5 = subplot('Position',[0.02 0.02 0.18 0.28]);
SP6 = subplot('Position',[0.22 0.02 0.21 0.28]);
SP7 = subplot('Position',[0.55 0.02 0.18 0.28]);
SP8 = subplot('Position',[0.75 0.02 0.21 0.28]);

CheckingFigure = figure('Name','CheckingFigure');
colormap(jet)
%% Panel A
StimInd = find(MasterMETA(PVExampleCell).StimFlag==1); StimInd = StimInd(PVStimUpStates);
StimUp = squeeze(MasterMETA(PVExampleCell).UpStates(1,StimInd,:));
PreUpState(1,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,StimInd(1)-1,12550:12950))';
PreUpState(2,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,StimInd(2)-1,12550:12950))';
PreUpState(3,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,StimInd(3)-1,12550:12950))';
ISI = 0;

figure(MainFigure)
subplot(SP2)
hold on
plot([PreUpState(1,:) StimUp(1,1:10500)],'Color',[0.1 0.1 0.1],'LineWidth',2.0)
plot([PreUpState(2,:) StimUp(2,1:10500)],'Color',[0.35 0.35 0.35],'LineWidth',2.0)
plot([PreUpState(3,:) StimUp(3,1:10500)],'Color',[0.6 0.6 0.6],'LineWidth',2.0)

%rectangle('Position',[600 -72 1000 60],'FaceColor',[0.5 0.7 1 0.5],'EdgeColor','none');
for k = 1:15
    rectangle('Position',[800+ISI -75 10 70],'FaceColor',[0.5 0.7 1 0.75],'EdgeColor','none');
    ISI = ISI+30;
end

line([4410 5410],[-62 -62],'Color','k','LineWidth',2.0)
line([5410 5410],[-62 -52],'Color','k','LineWidth',2.0)
text(4675,-63.5,'500ms','FontWeight','bold','FontSize',12)
text(5510,-53,'10mV','FontWeight','bold','FontSize',12,'Rotation',270)

hold off
xlim([0 8500])
ylim([-80 -30])
box off
axis off


%% Panel B
NonStimInd = find(MasterMETA(PVExampleCell).StimFlag==0); NonStimInd = NonStimInd(PVNonUpStates);
NonStimUp = squeeze(MasterMETA(PVExampleCell).UpStates(1,NonStimInd,:));
PreUpState(1,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,NonStimInd(1),12550:12950))';
PreUpState(2,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,NonStimInd(2),12550:12950))';
PreUpState(3,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,NonStimInd(3),12550:12950))';

figure(MainFigure)
subplot(SP1)
hold on
plot([PreUpState(1,:) NonStimUp(1,1:10500)],'Color',[0.1 0.1 0.1],'LineWidth',2.0)
plot([PreUpState(2,:) NonStimUp(2,1:10500)],'Color',[0.35 0.35 0.35],'LineWidth',2.0)
plot([PreUpState(3,:) NonStimUp(3,1:10500)],'Color',[0.6 0.6 0.6],'LineWidth',2.0)
rectangle('Position',[800 -75 450 40],'FaceColor','none','EdgeColor','k','LineStyle','--','LineWidth',2.0);
hold off
xlim([0 8500])
ylim([-80 -30])
box off
axis off

%% Panel C - PV positive control

load(PVPosCtrl,'CH1','CH3')
PVChetaTrace = CH1{3}(3,81368:81950);
ISI = 0;

subplot(SP3)
plot(PVChetaTrace,'r','LineWidth',2.0)
hold on
for k = 1:15
    rectangle('Position',[50+ISI -75 10 100],'FaceColor',[0.5 0.7 1 0.45],'EdgeColor','none');
    ISI = ISI+30;
end
line([380 580],[-76.5 -76.5],'Color','k','LineWidth',2.0)
line([580 580],[-76.5 -56.5],'Color','k','LineWidth',2.0)
text(400,-78,'100ms','FontWeight','bold','FontSize',12)
text(605,-62.5,'20mV','FontWeight','bold','FontSize',12,'Rotation',270)
hold off
ylim([-78 -38])
box off
axis off

%% Panel D - SST positive control
load(SSTPosCtrl,'CH1','CH3')
SSTChetaTrace = CH1{5}(2,4812:5912);
ISI = 0;

subplot(SP4)
plot(SSTChetaTrace,'b','LineWidth',2.0)
hold on
for k = 1:25
    rectangle('Position',[50+ISI -72 20 100],'FaceColor',[0.5 0.7 1 0.45],'EdgeColor','none');
    ISI = ISI+40;
end
line([910 1110],[-72.5 -72.5],'Color','k','LineWidth',2.0)
line([1110 1110],[-72.5 -52.5],'Color','k','LineWidth',2.0)
text(900,-74,'100ms','FontWeight','bold','FontSize',12)
text(1160,-58.5,'20mV','FontWeight','bold','FontSize',12,'Rotation',270)
hold off
ylim([-73 -32])
box off
axis off

%% Calculating Duration, Voltage Change, and Firing Rate for PV
for i = 1:length(MasterMETA) %Organazing the Durations of the UpStates
    FIRINGRATE = [];
    StimFRInd = [];
    
    PassNon = nan;
    PassStim = nan;
    if length(MasterMETA(i).UpDurPY)>9 && i~=10 && i~=11
        StimIndex       = MasterMETA(i).StimFlag==1;
        NonStimIndex    = MasterMETA(i).StimFlag==0;
        
        NonStimUpDur    = (MasterMETA(i).UpDurPY(1,NonStimIndex,:))*(Param.dt/1000);
        StimUpDur       = (MasterMETA(i).UpDurPY(1,StimIndex,:))*(Param.dt/1000);
        PassNon         = NonStimUpDur(NonStimUpDur>MinDuration);
        PassStim        = StimUpDur(StimUpDur>MinDuration);
        if length(PassNon)<5||length(PassStim)<5
            PassNon = nan;
            PassStim = nan;
        end
    end
    PVMETA(i).PassNonStimUpDur = PassNon;
    PVMETA(i).MedPassNonStimUpDur = median(PassNon);
    PVMETA(i).PassStimUpDur = PassStim;
    PVMETA(i).MedPassStimUpDur = median(PassStim);
    PVMETA(i).FailedStimUpDur  = median(PassStim(PassStim>FailMinLen));

    if length(MasterMETA(i).UpDurPY)>9 && i~=10 && i~=11 %Finding the Volage during stimmed and non stimmed
        CH1     = squeeze(MasterMETA(i).UpStates(1,:,:));
        CH3     = squeeze(MasterMETA(i).UpStates(3,:,:));
        UpDurPY = MasterMETA(i).UpDurPY;
        chunk   = zeros(size(CH3,1),(prestimsize+poststimsize));
        StimInd = zeros(size(CH3,1),1);

        
        %% Getting chunk centered around the first stimulus or 250 ms from onset (in non stimulated)
        for j = 1:size(CH3,1)
            if UpDurPY(j)>minlength
                ind = find(CH3(j,1:prestimsize*3)>10,1);
                if ~isempty(ind) && ind>500
                    chunk(j,:) = CH1(j,ind-prestimsize:ind+(poststimsize)-1);
                    StimInd(j) = ind;
                elseif isempty(ind)
                    chunk(j,:) = CH1(j,1+((Param.preEventTime+stimdelay)-prestimsize):(Param.preEventTime+stimdelay)+(poststimsize));
                end
            end
        end
        
%         offvoltwin       = sort([startlight+voltwinsize startlight]); 
        BadVarName       = [(voltstimlen/2)-(voltwinsize/2) (voltstimlen/2)+(voltwinsize/2)];
        offvoltwin       = BadVarName+startlight;
        filtchunk        = medfilt1(chunk',filtorder,'truncate')';
        
        AllStimVolt      = filtchunk(StimIndex&any(filtchunk,2)',:);
        StimVolt         = mean(AllStimVolt,1);
        AllNonVolt       = filtchunk(NonStimIndex&any(filtchunk,2)',:);
        NonVolt          = mean(AllNonVolt,1);
        
%         StimEndStim      = mean(StimVolt(offvoltwin(1)+voltstimlen:offvoltwin(2)+voltstimlen));
%         NonEndStim       = mean(NonVolt(offvoltwin(1)+voltstimlen:offvoltwin(2)+voltstimlen));
        StimEndStim      = mean(StimVolt(offvoltwin(1):offvoltwin(2)));
        NonEndStim       = mean(NonVolt(offvoltwin(1):offvoltwin(2)));
        StimPreStim      = mean(StimVolt(offvoltwin(1):offvoltwin(2)));
        NonPreStim       = mean(NonVolt(offvoltwin(1):offvoltwin(2)));
        
        FirstStimInd     = find(StimInd>0,1);
        FirstStimStart   = StimInd(FirstStimInd);
        stimchannel      = CH3(FirstStimInd,FirstStimStart-prestimsize:FirstStimStart+(poststimsize)-1);
        
        figure(CheckingFigure)
        plot(StimVolt,'b','LineWidth',2.0)
        hold on
        plot(NonVolt,'k','LineWidth',2.0)
        plot(AllStimVolt','Color',[0 0 0.5 0.5])
        plot(AllNonVolt','Color',[0.25 0.25 0.25 0.5])
        rectangle('Position',[offvoltwin(1) -85 (offvoltwin(2)-offvoltwin(1)) 85],'FaceColor',[0.5 0.7 1 0.5])
        rectangle('Position',[offvoltwin(1)+voltstimlen -85 (offvoltwin(2)-offvoltwin(1)) 85],'FaceColor',[0.5 0.7 1 0.5])
        plot(stimchannel(1,:),'LineWidth',2.5)
        hold off
        ylim([-85 10])
%         waitforbuttonpress
        clf
        
        if size(AllStimVolt,1)<5||size(AllNonVolt,1)<5
            StimPreStim    = nan;
            StimEndStim    = nan;
            NonPreStim     = nan;
            NonEndStim     = nan;
        end
        
        PVMETA(i).EndStim  = [StimEndStim; NonEndStim];
        PVMETA(i).PreStim  = [StimPreStim; NonPreStim];
        
 %% Getting Pyramidal Firing Rate     This one is big :/
        for j = 1:size(CH1,1)
            ind = find(CH3(j,1:prestimsize*3)>10,1);
            if UpDurPY(j)>minlengthFR && ((~isempty(ind) && ind>500) || isempty(ind))
                %Selecting Timeframes
                if ~isempty(ind) && ind>500
                    prewindow       = Param.preEventTime:ind;
                    swindow         = ind:ind+stimlen;
                    postwindow      = ind+stimlen:Param.preEventTime+UpDurPY(j);
                    
                    preS            = CH1(j,prewindow);
                    S               = CH1(j,swindow);
                    postS           = CH1(j,round(postwindow));                    
                    
                    if isempty(postwindow)|| length(postwindow)<3
                        postwindow = 0;
                        postS = [0 0 0];
                    end
                    
                    [~,preStimes]   = findpeaks(preS,'MinPeakHeight',spikethr);
                    [~,Stimes]      = findpeaks(S,'MinPeakHeight',spikethr);
                    [~,postStimes]  = findpeaks(postS,'MinPeakHeight',spikethr);
                    
                    preStimes       = preStimes +Param.preEventTime;
                    Stimes          = Stimes +ind;
                    postStimes      = postStimes +ind+stimlen;
                    
                    preSFR          = length(preStimes)/((length(prewindow))*(dt/1000));
                    SFR             = length(Stimes)/((length(swindow))*(dt/1000)); 
                    postSFR         = length(postStimes)/((length(postwindow))*(dt/1000));
                    
                elseif isempty(ind)
                    prewindow       = Param.preEventTime:Param.preEventTime+stimdelay;
                    swindow         = Param.preEventTime+stimdelay:Param.preEventTime+stimdelay+stimlen;
                    postwindow      = Param.preEventTime+stimdelay+stimlen:Param.preEventTime+UpDurPY(j);
                    
                    preS            = CH1(j,prewindow);
                    S               = CH1(j,swindow);
                    postS           = CH1(j,round(postwindow));
                    
                    if isempty(postwindow) || length(postwindow)<3
                        postwindow = 0;
                        postS = [0 0 0];
                    end
                    
                    [~,preStimes]   = findpeaks(preS,'MinPeakHeight',spikethr);
                    [~,Stimes]      = findpeaks(S,'MinPeakHeight',spikethr);
                    [~,postStimes]  = findpeaks(postS,'MinPeakHeight',spikethr);
                    
                    preStimes       = preStimes +Param.preEventTime;
                    Stimes          = Stimes +Param.preEventTime+stimdelay;
                    postStimes      = postStimes +Param.preEventTime+stimdelay+stimlen;
                    
                    preSFR          = length(preStimes)/((length(prewindow))*(dt/1000));
                    SFR             = length(Stimes)/((length(swindow))*(dt/1000)); 
                    postSFR         = length(postStimes)/((length(postwindow))*(dt/1000)); 
                    
                end                
                
                fr = [preSFR; SFR; postSFR];
                StimFRInd = [StimFRInd StimIndex(j)];
                FIRINGRATE = [FIRINGRATE fr];
                
%                 figure(CheckingFigure)
%                 plot(CH1(j,:))
%                 hold on
%                 rectangle('Position',[prewindow(1) -80 length(prewindow) 100],'FaceColor',[0 0 0 0.25],'EdgeColor','none')
%                 rectangle('Position',[swindow(1) -80 length(swindow) 100],'FaceColor',[0 1 0 0.25],'EdgeColor','none')
%                 rectangle('Position',[postwindow(1) -80 length(postwindow) 100],'FaceColor',[0 0 0 0.25],'EdgeColor','none')
%                 scatter(preStimes,repmat(spikethr,1,length(preStimes)),[],'w','filled')
%                 scatter(Stimes,repmat(spikethr,1,length(Stimes)),[],'y','filled')
%                 scatter(postStimes,repmat(spikethr,1,length(postStimes)),[],'w','filled')
%                 hold off
%                 ylim([-80 20])
%                 waitforbuttonpress
%                 clf
            end
            
        end        
        PVMETA(i).StimFR = FIRINGRATE(:,StimFRInd==1);
        PVMETA(i).MedStimFR = median(FIRINGRATE(:,StimFRInd==1),2);
        PVMETA(i).NonStimFR = FIRINGRATE(:,StimFRInd==0);
        PVMETA(i).MedNonStimFR = median(FIRINGRATE(:,StimFRInd==0),2);
        
    end
end

%% Panel E
PVNonDur  = [PVMETA(:).MedPassNonStimUpDur];
PVStimDur = [PVMETA(:).MedPassStimUpDur];

[p1,h1] = signrank(PVNonDur,PVStimDur);

figure(MainFigure)
subplot(SP5)
boxplot([PVNonDur' PVStimDur'],'Labels',{'Off','On'},'Symbol','o','Width',0.55)
% title('PV','FontSize',16,'FontWeight','bold')
ylim([0 5.5])
xlim([0.45 3.2])
A = gca;
set(A,'LineWidth',2,'FontSize',12,'FontWeight','bold')
set(A.Children.Children,'Color','k','LineWidth',2)
set(A.Children.Children(1:2),'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1,'MarkerSize',3)
set(A.Children.Children(11:14),'LineStyle','-')
color = {[1 0 0],'k'};
H = findobj(gca,'Tag','Box');
for j=1:length(H)
   patch(get(H(j),'XData'),get(H(j),'YData'),color{j},'FaceAlpha',.5);
end
hold on
line(repmat([1;2],1,size(PVNonDur,2)),[PVNonDur;PVStimDur],'Color',[0 0 0 0.25],'LineWidth',0.5)
hold off
ylabel('Pyr Up State Duration (s)','FontWeight','bold','FontSize',12)
box off
 
if h1==1
    hold on
    line([1 2],[5.1 5.1],'LineWidth',2,'Color','k')
    text(1.25,5.4,['\ast' '\ast'],'FontWeight','bold','FontSize',16)
end
%% Calculating Duration, Voltage Change, and Firing Rate for SOM
clear MasterMETA
cd('C:\Users\IzquierdoLab\Desktop\Juan Luis\PV_SST_Opto_Study_19\OptoUpStateData\SOM-Cheta\DoubleTrans')
load('SSTCheta2TransMETA.mat')
MasterMETA = MasterMETA(SSTGroup);
% voltstimlen = 490/dt;

for i = 1:length(MasterMETA)
    FIRINGRATE = [];
    StimFRInd = [];
   
    PassNon = nan;
    PassStim = nan;
    if length(MasterMETA(i).UpDurPY)>9
        StimIndex       = MasterMETA(i).StimFlag==1;
        NonStimIndex    = MasterMETA(i).StimFlag==0;
        
        NonStimUpDur    = (MasterMETA(i).UpDurPY(1,NonStimIndex,:))*(Param.dt/1000);
        StimUpDur       = (MasterMETA(i).UpDurPY(1,StimIndex,:))*(Param.dt/1000);
        PassNon         = NonStimUpDur(NonStimUpDur>MinDuration);
        PassStim        = StimUpDur(StimUpDur>MinDuration);
        if length(PassNon)<5||length(PassStim)<5
            PassNon = nan;
            PassStim = nan;
        end
    end
    SSTMETA(i).PassNonStimUpDur = PassNon;
    SSTMETA(i).MedPassNonStimUpDur = median(PassNon);
    SSTMETA(i).PassStimUpDur = PassStim;
    SSTMETA(i).MedPassStimUpDur = median(PassStim);
    SSTMETA(i).FailedStimUpDur  = median(PassStim(PassStim>FailMinLen));
    
    if length(MasterMETA(i).UpDurPY)>9  %Finding the Volage during stimmed and non stimmed       
        %% Getting chunk centered around the first stimulus or 250 ms from onset (in non stimulated)
        CH1     = squeeze(MasterMETA(i).UpStates(1,:,:));
        CH3     = squeeze(MasterMETA(i).UpStates(3,:,:));
        UpDurPY = MasterMETA(i).UpDurPY;
        chunk   = zeros(size(CH3,1),(prestimsize+poststimsize));
        StimInd = zeros(size(CH3,1),1);
                
        %% Getting chunk centered around the first stimulus or 250 ms from onset (in non stimulated)
        for j = 1:size(CH3,1)
            if UpDurPY(j)>minlength
                ind = find(CH3(j,1:prestimsize*3)>10,1);
                if ~isempty(ind) && ind>500
                    chunk(j,:) = CH1(j,ind-prestimsize:ind+(poststimsize)-1);
                    StimInd(j) = ind;
                elseif isempty(ind)
                    chunk(j,:) = CH1(j,1+((Param.preEventTime+stimdelay)-prestimsize):(Param.preEventTime+stimdelay)+(poststimsize));
                end
            end
        end
        
%         offvoltwin       = sort([startlight+voltwinsize startlight]); 
        BadVarName       = [(voltstimlen/2)-(voltwinsize/2) (voltstimlen/2)+(voltwinsize/2)];
        offvoltwin       = BadVarName+startlight;
        filtchunk        = medfilt1(chunk',filtorder,'truncate')';
        
        AllStimVolt      = filtchunk(StimIndex&any(filtchunk,2)',:);
        StimVolt         = mean(AllStimVolt,1);
        AllNonVolt       = filtchunk(NonStimIndex&any(filtchunk,2)',:);
        NonVolt          = mean(AllNonVolt,1);
        
%         StimEndStim      = mean(StimVolt(offvoltwin(1)+voltstimlen:offvoltwin(2)+voltstimlen));
%         NonEndStim       = mean(NonVolt(offvoltwin(1)+voltstimlen:offvoltwin(2)+voltstimlen));
        StimEndStim      = mean(StimVolt(offvoltwin(1):offvoltwin(2)));
        NonEndStim       = mean(NonVolt(offvoltwin(1):offvoltwin(2)));
        StimPreStim      = mean(StimVolt(offvoltwin(1):offvoltwin(2)));
        NonPreStim       = mean(NonVolt(offvoltwin(1):offvoltwin(2)));
        
        FirstStimInd     = find(StimInd>0,1);
        FirstStimStart   = StimInd(FirstStimInd);
        stimchannel      = CH3(FirstStimInd,FirstStimStart-prestimsize:FirstStimStart+(poststimsize)-1);
        
        figure(CheckingFigure)
        plot(StimVolt,'b','LineWidth',2.0)
        hold on
        plot(NonVolt,'k','LineWidth',2.0)
        plot(AllStimVolt','Color',[0 0 0.5 0.5])
        plot(AllNonVolt','Color',[0.25 0.25 0.25 0.5])
        rectangle('Position',[offvoltwin(1) -85 (offvoltwin(2)-offvoltwin(1)) 85],'FaceColor',[0.5 0.7 1 0.5])
        rectangle('Position',[offvoltwin(1)+voltstimlen -85 (offvoltwin(2)-offvoltwin(1)) 85],'FaceColor',[0.5 0.7 1 0.5])
        plot(stimchannel(1,:),'LineWidth',2.5)
        hold off
        ylim([-85 10])
%         waitforbuttonpress
        clf
        
        if size(AllStimVolt,1)<1||size(AllNonVolt,1)<1
            StimPreStim    = nan;
            StimEndStim    = nan;
            NonPreStim     = nan;
            NonEndStim     = nan;
        end
        
        SSTMETA(i).EndStim  = [StimEndStim; NonEndStim];
        SSTMETA(i).PreStim  = [StimPreStim; NonPreStim];
        %% Getting Pyramidal Firing Rate     This one is big :/
        for j = 1:size(CH1,1)
            ind = find(CH3(j,1:prestimsize*3)>10,1);
            if UpDurPY(j)>minlengthFR && ((~isempty(ind) && ind>500) || isempty(ind))
                %Selecting Timeframes
                if ~isempty(ind) && ind>500
                    prewindow       = Param.preEventTime:ind;
                    swindow         = ind:ind+stimlen;
                    postwindow      = ind+stimlen:Param.preEventTime+UpDurPY(j);
                    
                    preS            = CH1(j,prewindow);
                    S               = CH1(j,swindow);
                    postS           = CH1(j,round(postwindow));                    
                    
                    if isempty(postwindow)|| length(postwindow)<3
                        postwindow = 0;
                        postS = [0 0 0];
                    end
                    
                    [~,preStimes]   = findpeaks(preS,'MinPeakHeight',spikethr);
                    [~,Stimes]      = findpeaks(S,'MinPeakHeight',spikethr);
                    [~,postStimes]  = findpeaks(postS,'MinPeakHeight',spikethr);
                    
                    preStimes       = preStimes +Param.preEventTime;
                    Stimes          = Stimes +ind;
                    postStimes      = postStimes +ind+stimlen;
                    
                    preSFR          = length(preStimes)/((length(prewindow))*(dt/1000));
                    SFR             = length(Stimes)/((length(swindow))*(dt/1000)); 
                    postSFR         = length(postStimes)/((length(postwindow))*(dt/1000));
                    
                elseif isempty(ind)
                    prewindow       = Param.preEventTime:Param.preEventTime+stimdelay;
                    swindow         = Param.preEventTime+stimdelay:Param.preEventTime+stimdelay+stimlen;
                    postwindow      = Param.preEventTime+stimdelay+stimlen:Param.preEventTime+UpDurPY(j);
                    
                    preS            = CH1(j,prewindow);
                    S               = CH1(j,swindow);
                    postS           = CH1(j,round(postwindow));
                    
                    if isempty(postwindow)|| length(postwindow)<3
                        postwindow = 0;
                        postS = [0 0 0];
                    end
                    
                    [~,preStimes]   = findpeaks(preS,'MinPeakHeight',spikethr);
                    [~,Stimes]      = findpeaks(S,'MinPeakHeight',spikethr);
                    [~,postStimes]  = findpeaks(postS,'MinPeakHeight',spikethr);
                    
                    preStimes       = preStimes +Param.preEventTime;
                    Stimes          = Stimes +Param.preEventTime+stimdelay;
                    postStimes      = postStimes +Param.preEventTime+stimdelay+stimlen;
                    
                    preSFR          = length(preStimes)/((length(prewindow))*(dt/1000));
                    SFR             = length(Stimes)/((length(swindow))*(dt/1000)); 
                    postSFR         = length(postStimes)/((length(postwindow))*(dt/1000)); 
                    
                end                
                
                fr = [preSFR; SFR; postSFR];
                StimFRInd = [StimFRInd StimIndex(j)];
                FIRINGRATE = [FIRINGRATE fr];
                
%                 figure(CheckingFigure)
%                 plot(CH1(j,:))
%                 hold on
%                 rectangle('Position',[prewindow(1) -80 length(prewindow) 100],'FaceColor',[0 0 0 0.25],'EdgeColor','none')
%                 rectangle('Position',[swindow(1) -80 length(swindow) 100],'FaceColor',[0 1 0 0.25],'EdgeColor','none')
%                 rectangle('Position',[postwindow(1) -80 length(postwindow) 100],'FaceColor',[0 0 0 0.25],'EdgeColor','none')
%                 scatter(preStimes,repmat(spikethr,1,length(preStimes)),[],'w','filled')
%                 scatter(Stimes,repmat(spikethr,1,length(Stimes)),[],'y','filled')
%                 scatter(postStimes,repmat(spikethr,1,length(postStimes)),[],'w','filled')
%                 hold off
%                 ylim([-80 20])
%                 waitforbuttonpress
%                 clf
            end
            
        end        
        SSTMETA(i).StimFR = FIRINGRATE(:,StimFRInd==1);
        SSTMETA(i).MedStimFR = median(FIRINGRATE(:,StimFRInd==1),2);
        SSTMETA(i).NonStimFR = FIRINGRATE(:,StimFRInd==0);
        SSTMETA(i).MedNonStimFR = median(FIRINGRATE(:,StimFRInd==0),2);
    end
end

%% Panel F
SSTNonDur = [SSTMETA(:).MedPassNonStimUpDur];
SSTStimDur = [SSTMETA(:).MedPassStimUpDur];

[p2,h2] = signrank(SSTNonDur,SSTStimDur);

figure(MainFigure)
subplot(SP7)
boxplot([SSTNonDur' SSTStimDur'],'Labels',{'Off','On'},'Symbol','o','Width',0.55)
% title('SOM','FontSize',16,'FontWeight','bold')
ylim([0 5.5])
xlim([0.45 3.2])
A = gca;
set(A,'LineWidth',2,'FontSize',12,'FontWeight','bold')
set(A.Children.Children,'Color','k','LineWidth',2)
set(A.Children.Children(1:2),'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1,'MarkerSize',3)
set(A.Children.Children(11:14),'LineStyle','-')
color = {[0 0 1],'k'};
H = findobj(gca,'Tag','Box');
for j=1:length(H)
   patch(get(H(j),'XData'),get(H(j),'YData'),color{j},'FaceAlpha',.5);
end
hold on
line(repmat([1;2],1,size(SSTNonDur,2)),[SSTNonDur;SSTStimDur],'Color',[0 0 0 0.25],'LineWidth',0.5)
hold off
ylabel('Pyr Up State Duration (s)','FontWeight','bold','FontSize',12)
box off

if h2==1
    hold on
    line([1 2],[5 5],'LineWidth',2,'Color','k')
    text(1.4,5.3,['\ast' '\ast'],'FontWeight','bold','FontSize',16)
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Change in Voltage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Graphics
% PanelC = figure('Name','MainFigure');
% set(gcf,'color','w');
% ScreenSize = get(groot,'ScreenSize');
% PanelC.Position = [50 100 ScreenSize(3)/1.8 ScreenSize(4)/1.3]; %This is to Force the figure to have a vertical orientation
% 
% SP1 = subplot('Position',[0.05 0.05 0.4 0.3]);
% SP2 = subplot('Position',[0.55 0.05 0.4 0.3]);


%% Getting Data
EndStimPV   = [PVMETA(:).EndStim]; %The Stim is row 1 and Non-Stim is row 2
PreStimPV   = [PVMETA(:).PreStim];
EndStimSST  = [SSTMETA(:).EndStim];
PreStimSST  = [SSTMETA(:).PreStim];

subplot(SP6)
boxplot([EndStimPV(2,:)' EndStimPV(1,:)'],'Labels',{'Off','On'},'Symbol','o','Width',0.55)
% title('Difference in Voltage','FontSize',12)
ylabel('Voltage (mV)')
A = gca;
set(A,'LineWidth',2,'FontSize',12,'FontWeight','bold','YAxisLocation','right')
set(A.Children.Children,'Color','k','LineWidth',2)
set(A.Children.Children(1:2),'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1,'MarkerSize',3)
set(A.Children.Children(11:14),'LineStyle','-')
color = {[1 0 0],'k'};
H = findobj(gca,'Tag','Box');
for j=1:length(H)
   patch(get(H(j),'XData'),get(H(j),'YData'),color{j},'FaceAlpha',.5);
end
hold on
line(repmat([1;2],1,size(EndStimPV,2)),flipud(EndStimPV),'Color',[0 0 0 0.25],'LineWidth',0.5)
hold off
ylim([-70 -40])
xlim([0 2.7])
box off

[p33,h33] = signrank(EndStimPV(2,:),EndStimPV(1,:));

if h33 == 1
    hold on
    line([1 2],[-42 -42],'Color','k','LineWidth',2) 
    text(1.21,-40.6,'\ast \ast','FontWeight','bold','FontSize',16)
end

subplot(SP8)
boxplot([EndStimSST(2,:)' EndStimSST(1,:)'],'Labels',{'Off','On'},'Symbol','o','Width',0.55)
% title('Difference in Voltage','FontSize',12)
ylabel('Voltage (mV)')
A = gca;
set(A,'LineWidth',2,'FontSize',12,'FontWeight','bold','YAxisLocation','right')
set(A.Children.Children,'Color','k','LineWidth',2)
set(A.Children.Children(1:2),'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1,'MarkerSize',3)
set(A.Children.Children(11:14),'LineStyle','-')
color = {[0 0 1],'k'};
H = findobj(gca,'Tag','Box');
for j=1:length(H)
   patch(get(H(j),'XData'),get(H(j),'YData'),color{j},'FaceAlpha',.5);
end
hold on
line(repmat([1;2],1,size(EndStimSST,2)),flipud(EndStimSST),'Color',[0 0 0 0.25],'LineWidth',0.5)
hold off
ylim([-70 -40])
xlim([-0.1 2.7])
box off

[p44,h44] = signrank(EndStimSST(2,:),EndStimSST(1,:));

if h44 == 1
    hold on
    line([1 2],[-42 -42],'Color','k','LineWidth',2) 
    text(1.21,-40.6,'\ast \ast','FontWeight','bold','FontSize',14)
end

cd('C:\Users\IzquierdoLab\Desktop\Juan Luis\PV_SST_Opto_Study_19\FigureRepo')
print '-dtiffn' OptoFigureCheta_WithSpikes

% % FiringRate Section
% PVStimFR     = [PVMETA(:).MedStimFR];
% PVNonStimFR  = [PVMETA(:).MedNonStimFR];
% 
% SOMStimFR    = [SSTMETA(:).MedStimFR];
% SOMNonStimFR = [SSTMETA(:).MedNonStimFR];
% 
% [p8,h8] = signrank(PVStimFR(2,:),PVNonStimFR(2,:));
% [p9,h9] = signrank(SOMStimFR(2,:),SOMNonStimFR(2,:));
% 
% figure(CheckingFigure)
% subplot(1,2,1)
% boxplot([PVStimFR(2,:)' PVNonStimFR(2,:)'],'Labels',{'Stim','NonStim'})
% title('PV PYFR')
% subplot(1,2,2)
% boxplot([SOMStimFR(2,:)' SOMNonStimFR(2,:)'],'Labels',{'Stim','NonStim'})
% title('SOM PYFR')
% 
% % Robust Stimmed UpState duration
% PVRobustDur  = [PVMETA(:).FailedStimUpDur];
% SOMRobustDur = [SSTMETA(:).FailedStimUpDur];
% 
% [p11,h11] = signrank(PVNonDur,PVRobustDur);
% [p22,h22] = signrank(SSTNonDur,SOMRobustDur);
% 
% figure
% subplot(1,2,1)
% boxplot([PVNonDur' PVRobustDur'],'Labels',{'NonStim','Robust'})
% title('PV DUR')
% subplot(1,2,2)
% boxplot([SSTNonDur' SOMRobustDur'],'Labels',{'NonStim','Robust'})
% title('SOM DUR')
% 
% % Trash
% ALLNSTIMUPDUR = [];
% for i = 1:length(PVMETA)
% AllNStim = PVMETA(i).PassStimUpDur/max(PVMETA(i).PassStimUpDur);
% 
% ALLNSTIMUPDUR = [ALLNSTIMUPDUR AllNStim];
% end

