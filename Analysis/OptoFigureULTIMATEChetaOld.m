close all
clear
warning('off','signal:findpeaks:largeMinPeakHeight')

cd('C:\Users\BuonoLab\Desktop\Juan Luis\MainStudy\OptogeneticFigureforGRant')
load('PVChetaMETA.mat')

%Variables for the Duration of UpState
PVExampleCell = 4;
PVStimUpStates = [17,21,28];
PVNonUpStates = [34,37,50];
SOMExample = 0;
MinDuration = 0.25;
FailMinLen = 0.750;

%Variables for the Difference in Voltage
dt = Param.dt;
prestimsize = 275/dt;
poststimsize = 600/dt;
stimlen = 500/dt;
startlight = prestimsize;
minlength = 250/dt;
stimdelay = 250/dt;
voltwinsize = -100/dt;
filtorder = 5/dt;

%Variables for Firing Rate
spikethr = -25;
minlengthFR = 250/dt;


PVMETA = MasterMETA;

%% INIT GRAPHICS
MainFigure = figure('Name','MainFigure');
set(gcf,'color','w');
ScreenSize = get(groot,'ScreenSize');
MainFigure.Position = [50 30 ScreenSize(3)/1.6 ScreenSize(4)/1.1]; %This is to Force the figure to have a vertical orientation

SP1 = subplot('Position',[0.01 0.82 0.95 0.18]);
SP2 = subplot('Position',[0.01 0.64 0.95 0.18]);
SP3 = subplot('Position',[0.08 0.35 0.3 0.25]);
SP4 = subplot('Position',[0.58 0.35 0.3 0.25]);
SP5 = subplot('Position',[0.05 0.02 0.4 0.28]);
SP6 = subplot('Position',[0.55 0.02 0.4 0.28]);

CheckingFigure = figure('Name','CheckingFigure');

colormap(jet)
%% Panel A
StimInd = find(MasterMETA(PVExampleCell).StimFlag==1); StimInd = StimInd(PVStimUpStates);
StimUp = squeeze(MasterMETA(PVExampleCell).UpStates(1,StimInd,:));
PreUpState(1,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,StimInd(1)-1,12750:13650))';
PreUpState(2,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,StimInd(1)-1,7270:8170))';
PreUpState(3,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,StimInd(1)-1,6530:7430))';
ISI = 0;

figure(MainFigure)
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
xlim([0 6415])
ylim([-72 -43])
box off
axis off


%% Panel B
NonStimInd = find(MasterMETA(PVExampleCell).StimFlag==0); NonStimInd = NonStimInd(PVNonUpStates);
NonStimUp = squeeze(MasterMETA(PVExampleCell).UpStates(1,NonStimInd,:));
PreUpState(1,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,NonStimInd(1)-1,2800:3700))';
PreUpState(2,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,NonStimInd(2)-1,950:1850))';
PreUpState(3,:) = squeeze(MasterMETA(PVExampleCell).UpStates(1,NonStimInd(2)-1,6450:7350))';

figure(MainFigure)
subplot(SP1)
hold on
plot([PreUpState(1,:) NonStimUp(1,1:5500)],'Color',[0.1 0.1 0.1],'LineWidth',2.0)
plot([PreUpState(2,:) NonStimUp(2,1:5500)],'Color',[0.35 0.35 0.35],'LineWidth',2.0)
plot([PreUpState(3,:) NonStimUp(3,1:5500)],'Color',[0.6 0.6 0.6],'LineWidth',2.0)
rectangle('Position',[1500 -70 1000 27],'FaceColor','none','EdgeColor','k','LineStyle','--','LineWidth',2.0);
line([5410 6410],[-69 -69],'Color','k','LineWidth',2.0)
line([6410 6410],[-69 -59],'Color','k','LineWidth',2.0)
text(5675,-70.3,'500ms','FontWeight','bold','FontSize',12)
text(6510,-60,'10mV','FontWeight','bold','FontSize',12,'Rotation',270)
hold off
xlim([0 6415])
ylim([-71 -42.5])
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
               
        offvoltwin       = sort([startlight+voltwinsize startlight]);
        
        filtchunk        = medfilt1(chunk',filtorder,'truncate')';
        OffVolt          = mean(filtchunk(:,offvoltwin(1):offvoltwin(2)),2);
        OnVolt           = mean(filtchunk(:,offvoltwin(1)+stimlen:offvoltwin(2)+stimlen),2);
        
        StimOffVolt      = OffVolt(StimIndex&(any(OffVolt,2)'));
        NonOffVolt       = OffVolt(NonStimIndex&(any(OffVolt,2)'));
        StimOnVolt       = OnVolt(StimIndex&(any(OnVolt,2)'));
        NonOnVolt        = OnVolt(NonStimIndex&(any(OnVolt,2)'));
        
        FirstStimInd     = find(StimInd>0,1);
        FirstStimStart   = StimInd(FirstStimInd);
        stimchannel      = CH3(FirstStimInd,FirstStimStart-prestimsize:FirstStimStart+(poststimsize)-1);
        
        figure(CheckingFigure)
        plot(filtchunk(any(filtchunk,2),:)')
        hold on
        rectangle('Position',[offvoltwin(1) -85 (offvoltwin(2)-offvoltwin(1)) 85],'FaceColor',[0.5 0.7 1 0.5])
        rectangle('Position',[offvoltwin(1)+stimlen -85 (offvoltwin(2)-offvoltwin(1)) 85],'FaceColor',[0.5 0.7 1 0.5])
        plot(stimchannel(1,:),'LineWidth',2.5)
        hold off
        ylim([-85 10])
%         waitforbuttonpress
        clf
        
        if size(StimOffVolt,1)<5||size(NonOffVolt,1)<5
            StimOnVolt    = nan;
            StimOffVolt   = nan;
            NonOnVolt     = nan;
            NonOffVolt    = nan;
        end
        
        PVMETA(i).OffVolt  = [mean(StimOffVolt); mean(NonOffVolt)];
        PVMETA(i).OnVolt   = [mean(StimOnVolt); mean(NonOnVolt)];
        
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

%% Panel C
PVNonDur  = [PVMETA(:).MedPassNonStimUpDur];
PVStimDur = [PVMETA(:).MedPassStimUpDur];

[p1,h1] = signrank(PVNonDur,PVStimDur);

figure(MainFigure)
subplot(SP3)
boxplot([PVNonDur' PVStimDur'],'Labels',{'Off','On'},'Symbol','o')
title('PV','FontSize',16,'FontWeight','bold')
ylim([0 5.5])
A = gca;
set(A,'LineWidth',2,'FontSize',12,'FontWeight','bold')
set(A.Children.Children,'Color','k','LineWidth',2)
set(A.Children.Children(1:3),'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1,'MarkerSize',3)
color = {[1 0 0],'k'};
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
%% Calculating Duration, Voltage Change, and Firing Rate for SOM
clear MasterMETA
load('SOMCheta2TransMETA.mat')
SOMMETA = MasterMETA;

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
    SOMMETA(i).PassNonStimUpDur = PassNon;
    SOMMETA(i).MedPassNonStimUpDur = median(PassNon);
    SOMMETA(i).PassStimUpDur = PassStim;
    SOMMETA(i).MedPassStimUpDur = median(PassStim);
    SOMMETA(i).FailedStimUpDur  = median(PassStim(PassStim>FailMinLen));
    
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
        
        offvoltwin       = sort([startlight+voltwinsize startlight]);
        
        filtchunk        = medfilt1(chunk',filtorder,'truncate')';
        OffVolt          = mean(filtchunk(:,offvoltwin(1):offvoltwin(2)),2);
        OnVolt           = mean(filtchunk(:,offvoltwin(1)+stimlen:offvoltwin(2)+stimlen),2);
        
        StimOffVolt      = OffVolt(StimIndex&(any(OffVolt,2)'));
        NonOffVolt       = OffVolt(NonStimIndex&(any(OffVolt,2)'));
        StimOnVolt       = OnVolt(StimIndex&(any(OnVolt,2)'));
        NonOnVolt        = OnVolt(NonStimIndex&(any(OnVolt,2)'));
        
        FirstStimInd     = find(StimInd>0,1);
        FirstStimStart   = StimInd(FirstStimInd);
        stimchannel      = CH3(FirstStimInd,FirstStimStart-prestimsize:FirstStimStart+(poststimsize)-1);
        
        figure(CheckingFigure)
        plot(filtchunk(any(filtchunk,2),:)')
        hold on
        rectangle('Position',[offvoltwin(1) -85 (offvoltwin(2)-offvoltwin(1)) 85],'FaceColor',[0.5 0.7 1 0.5])
        rectangle('Position',[offvoltwin(1)+stimlen -85 (offvoltwin(2)-offvoltwin(1)) 85],'FaceColor',[0.5 0.7 1 0.5])
        plot(stimchannel(1,:),'LineWidth',2.5)
        hold off
        ylim([-85 10])
%         waitforbuttonpress
        clf
        
        if size(StimOffVolt,1)<5||size(NonOffVolt,1)<5
            StimOnVolt    = nan;
            StimOffVolt   = nan;
            NonOnVolt     = nan;
            NonOffVolt    = nan;
        end
        SOMMETA(i).OnVolt  = [mean(StimOffVolt); mean(NonOffVolt)];
        SOMMETA(i).OffVolt = [mean(StimOnVolt); mean(NonOnVolt)];
        
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
        SOMMETA(i).StimFR = FIRINGRATE(:,StimFRInd==1);
        SOMMETA(i).MedStimFR = median(FIRINGRATE(:,StimFRInd==1),2);
        SOMMETA(i).NonStimFR = FIRINGRATE(:,StimFRInd==0);
        SOMMETA(i).MedNonStimFR = median(FIRINGRATE(:,StimFRInd==0),2);
    end
end

%% Panel D
SOMNonDur = [SOMMETA(:).MedPassNonStimUpDur];
PVStimDur = [SOMMETA(:).MedPassStimUpDur];

[p2,h2] = signrank(SOMNonDur,PVStimDur);

figure(MainFigure)
subplot(SP4)
boxplot([SOMNonDur' PVStimDur'],'Labels',{'Off','On'},'Symbol','o')
title('SOM','FontSize',16,'FontWeight','bold')
ylim([0 5.5])
A = gca;
set(A,'LineWidth',2,'FontSize',12,'FontWeight','bold')
set(A.Children.Children,'Color','k','LineWidth',2)
set(A.Children.Children(1:4),'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1,'MarkerSize',3)
color = {[0 0 1],'k'};
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
OnPVVoltage   = [PVMETA(:).OnVolt];
OffPVVoltage  = [PVMETA(:).OffVolt];
OnSOMVoltage  = [SOMMETA(:).OnVolt];
OffSOMVoltage = [SOMMETA(:).OffVolt];

subplot(SP5)
boxplot([OffPVVoltage(2,:)' OffPVVoltage(1,:)' OnPVVoltage(2,:)' OnPVVoltage(1,:)'],'Labels',{'PreStim','PostStim','',''},'Symbol','o')
% title('Difference in Voltage','FontSize',12)
ylabel('Voltage (mV)')
A = gca;
set(A,'LineWidth',2,'FontSize',12,'FontWeight','bold','XTick',[1.5 3.5])
set(A.Children.Children,'Color','k','LineWidth',2)
set(A.Children.Children(1:4),'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1,'MarkerSize',3)
color = {[1 0 0],'k',[1 0 0],'k'};
H = findobj(gca,'Tag','Box');
for j=1:length(H)
   patch(get(H(j),'XData'),get(H(j),'YData'),color{j},'FaceAlpha',.5);
end
ylim([-71 -38])
box off

[p3,h3]   = signrank(OffPVVoltage(2,:),OffPVVoltage(1,:));
[p33,h33] = signrank(OnPVVoltage(2,:),OnPVVoltage(1,:));

if h33 == 1
    hold on
    line([3 4],[-45 -45],'Color','k','LineWidth',2) 
    text(3.32,-44,'\ast \ast','FontWeight','bold','FontSize',14)
end

subplot(SP6)
boxplot([OnSOMVoltage(2,:)' OnSOMVoltage(1,:)' OffSOMVoltage(2,:)' OffSOMVoltage(1,:)'],'Labels',{'PreStim','PostStim','',''},'Symbol','o')
% title('Difference in Voltage','FontSize',12)
ylabel('Voltage (mV)')
A = gca;
set(A,'LineWidth',2,'FontSize',12,'FontWeight','bold','XTick',[1.5 3.5])
set(A.Children.Children,'Color','k','LineWidth',2)
set(A.Children.Children(1:4),'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1,'MarkerSize',3)
color = {[0 0 1],'k',[0 0 1],'k'};
H = findobj(gca,'Tag','Box');
for j=1:length(H)
   patch(get(H(j),'XData'),get(H(j),'YData'),color{j},'FaceAlpha',.5);
end
ylim([-71 -38])
box off

[p4,h4]   = signrank(OffSOMVoltage(2,:),OffSOMVoltage(1,:));
[p44,h44] = signrank(OnSOMVoltage(2,:),OnSOMVoltage(1,:));

if h44 == 1
    hold on
    line([3 4],[-45 -45],'Color','k','LineWidth',2) 
    text(3.2,-43,'\ast \ast','FontWeight','bold','FontSize',14)
end

print '-dtiffn' OptoFigureCheta

%% FiringRate Section
PVStimFR     = [PVMETA(:).MedStimFR];
PVNonStimFR  = [PVMETA(:).MedNonStimFR];

SOMStimFR    = [SOMMETA(:).MedStimFR];
SOMNonStimFR = [SOMMETA(:).MedNonStimFR];

[p8,h8] = signrank(PVStimFR(2,:),PVNonStimFR(2,:));
[p9,h9] = signrank(SOMStimFR(2,:),SOMNonStimFR(2,:));

figure(CheckingFigure)
subplot(1,2,1)
boxplot([PVStimFR(2,:)' PVNonStimFR(2,:)'],'Labels',{'Stim','NonStim'})
title('PV PYFR')
subplot(1,2,2)
boxplot([SOMStimFR(2,:)' SOMNonStimFR(2,:)'],'Labels',{'Stim','NonStim'})
title('SOM PYFR')

%% Robust Stimmed UpState duration
PVRobustDur  = [PVMETA(:).FailedStimUpDur];
SOMRobustDur = [SOMMETA(:).FailedStimUpDur];

[p11,h11] = signrank(PVNonDur,PVRobustDur);
[p22,h22] = signrank(SOMNonDur,SOMRobustDur);

figure
subplot(1,2,1)
boxplot([PVNonDur' PVRobustDur'],'Labels',{'NonStim','Robust'})
title('PV DUR')
subplot(1,2,2)
boxplot([SOMNonDur' SOMRobustDur'],'Labels',{'NonStim','Robust'})
title('SOM DUR')

%% Trash
ALLNSTIMUPDUR = [];
for i = 1:length(PVMETA)
AllNStim = PVMETA(i).PassStimUpDur/max(PVMETA(i).PassStimUpDur);

ALLNSTIMUPDUR = [ALLNSTIMUPDUR AllNStim];
end

