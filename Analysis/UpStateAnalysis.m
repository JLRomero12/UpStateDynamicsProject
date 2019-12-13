tic
clear META
close all
warning('off','signal:findpeaks:largeMinPeakHeight')

graphics = 1; %graphics on and off switch

%Variables for finding the Upstates
thr = 5;             %Threshold from the resting of the cell (median of data)
SegLen = 10000;      %This is the maximun length that we expect the upstate to be. (25000 = 5s)
MaxIntEvent = 100;   %Max difference between the ThrCross and still be considered the same even (100 = 50ms)
MaxIntUpstate = 100; %Max difference between the ThrCross and still be considered the same even (100 = 50ms)
MinAboveTime = 1000; %Time that the event has to be above the threshold (Within the event) to count as Upstate (1000 = 500ms)
MinOnOffInt = 1000;  %Time that the event has to be above the threshold (from beginning to end) to count as Upstate (2500 = 500ms)
preEventTime = 200;  %Time given before the start of the upstate

VStimLag = 500;
VStimLen = 380;

Files = dir('S*.mat');
numFiles = length(Files);

META(1:numFiles) = struct();

if graphics
    % UpStateFigure = figure ('Name','UpStateFigure');
    DifferenceFigure = figure('Name','DifferenceFigure');
    AllTrialsFigure = figure('Name','AllTrialsFigure');
    FiringRateFigure = figure('Name','FiringRateFigure');
    XCorrFigure = figure('Name','XCorrFigure');
end

for file = 1:numFiles
    clearvars -EXCEPT META Files numFiles file thr SegLen MaxIntEvent MaxIntUpstate MinAboveTime MinOnOffInt preEventTime UpStateFigure ...
        DifferenceFigure AllTrialsFigure FiringRateFigure XCorrFigure graphics VStimLag VStimLen
    
    load(Files(file).name)
    
    UpState1 = [];
    UpState2 = [];
    UpState3 = [];
    UPDUR1 = [];
    STIMFLAG = [];
    STARTDIFF = [];
    UPSTATECORR = [];
    UPSTATECORRWINDOW = [];
    FIRINGRATE = [];
    STIMFIRINGRATE = [];
    VIRTSTIMFIRINGRATE = [];
    PEAKSHIFTS = [];
    STIMTIMES = [];
    AREA50 = [];
    AREA100 = [];
    
    if ~isempty(find(strcmp(CondOrder,'NEW_BMI_LIGHT_FF'),1))
        UpStateIndex = find(strcmp(CondOrder,'NEW_BMI_LIGHT_FF'));
    end
    
    for cond = UpStateIndex
        
        ch1 = CH1{cond};
        ch2 = CH2{cond};
        ch3 = CH3{cond};
        dt = round(SampleRate{cond},1);
        
        if graphics
            figure(AllTrialsFigure)
            if isempty(AllTrialsFigure.Children)
                imagesc(ch1);
            elseif ~isempty(AllTrialsFigure.Children)
                AllTrialsFigure.Children.Children.CData = ch1;
            end
        end
        
        for i = 1:size(ch1,1)
            %Find the Upstate
            tracePY = ch1(i,:);
            tracePV = ch2(i,:);
            %             figure(UpStateFigure)
            [~, ~, ~, ~, ~, upstate1, updur1, upstart1, upend1]...
                =MAKESEGMENTS(tracePY,thr,MaxIntEvent,MaxIntUpstate,SegLen,MinAboveTime,MinOnOffInt,preEventTime,0);
            
            UpState1 = [UpState1; upstate1];
            UPDUR1   = [UPDUR1 updur1];
            
            thrPY = mode(tracePY)+thr;
            thrPV = mode(tracePV)+thr;
            
            for j = 1:length(upstart1)
                window = upstart1(j)-preEventTime:upstart1(j)+preEventTime;
                MaxDiff = MaxIntEvent; %30 ms
                
                %Find the difference in the upstate time
                if ~isempty(find(tracePV(window)>(mode(tracePV)+thr),1)) && ~isempty(find(tracePY(window)>(mode(tracePY)+thr),1))...
                        && (mode(tracePV(window))<thrPV)
                    
                    AbovePV = find(tracePV(window)>(thrPV));
                    OnsetStartPV = (find(diff([-9999 AbovePV])>MaxDiff));
                    OffsetStartPV=AbovePV([OnsetStartPV(2:length(OnsetStartPV))-1 length(AbovePV)]);
                    OnsetStartPV = AbovePV(OnsetStartPV);
                    [~,n] = max(OffsetStartPV-OnsetStartPV);
                    UpStartPV = OnsetStartPV(n);
                    clear n
                    
                    AbovePY = find(tracePY(window)>(thrPY));
                    OnsetStartPY = (find(diff([-9999 AbovePY])>MaxDiff));
                    OffsetStartPY=AbovePY([OnsetStartPY(2:length(OnsetStartPY))-1 length(AbovePY)]);
                    OnsetStartPY = AbovePY(OnsetStartPY);
                    [~,n] = max(OffsetStartPY-OnsetStartPY);
                    UpStartPY = OnsetStartPY(n);
                    
                    startdiff = UpStartPY-UpStartPV;
                    STARTDIFF = [STARTDIFF startdiff];
                    
                    %Visualize Difference to Check
                    if graphics
                        figure(DifferenceFigure)
                        if isempty(DifferenceFigure.Children)
                            plot(tracePY(window),'LineWidth',2,'Color','k');
                            hold on
                            line([UpStartPY UpStartPY],[-70 10],'Color',[0 0 0 0.75],'LineStyle','--') %In the window of data made upstart1 is 501
                            plot(tracePV(window),'LineWidth',2,'Color','r')
                            line([UpStartPV UpStartPV],[-70 10],'Color',[1 0 0 0.75],'LineStyle','--')
                            line([1 length(window)],[thrPY thrPY],'LineStyle','--','Color','k')
                            line([1 length(window)],[thrPV thrPV],'LineStyle','--','Color','r')
                            hold off
                        elseif ~isempty(get(gca,'Children'))
                            DifferenceFigure.Children.Children(6).YData = tracePY(window);
                            DifferenceFigure.Children.Children(5).XData = [UpStartPY UpStartPY];
                            DifferenceFigure.Children.Children(4).YData = tracePV(window);
                            DifferenceFigure.Children.Children(3).XData = [UpStartPV UpStartPV];
                            DifferenceFigure.Children.Children(2).YData = [thrPY thrPY];
                            DifferenceFigure.Children.Children(1).YData = [thrPV thrPV];
                        end
                    end
                end
                
                %Making Vectors for Upstates in Ch2 and Ch3
                upstate2(j,:) = ch2(i,upstart1(j)-preEventTime:upstart1(j)+SegLen);
                UpState2 = [UpState2; upstate2(j,:)];
                
                upstate3 = ch3(i,upstart1(j)-preEventTime:upstart1(j)+SegLen);
                UpState3 = [UpState3; upstate3];
                
                %Finding Stimulation
                stimtimes = find(upstate3>10);
                if ~isempty(stimtimes)
                    if stimtimes(1)>preEventTime && stimtimes(1)<(preEventTime+updur1(j))
                        stimtimes = stimtimes(stimtimes<(preEventTime+updur1(j))); 
                        stimtimes = [stimtimes(1);stimtimes(end)];
                    elseif stimtimes (1)<preEventTime || stimtimes(1)>(preEventTime+updur1(j))
                        stimtimes = [];
                    end
                end
                stimflag = ~isempty(stimtimes);
                STIMTIMES = [STIMTIMES stimtimes];
                STIMFLAG = [STIMFLAG stimflag];
                
                %Correlation of Filtered PY Upstate and PV Upstate
                FiltUpStatePY = medfilt1(ch1(i,upstart1(j)-preEventTime:upend1(j)+preEventTime),50,'truncate');
                FiltUpStatePV = medfilt1(ch2(i,upstart1(j)-preEventTime:upend1(j)+preEventTime),50,'truncate');
                c = corr(FiltUpStatePY',FiltUpStatePV');
                UPSTATECORR = [UPSTATECORR c];
                
                %Correlation of Filtered PY Upstate and PV Upstate for full
                %window (with resting)
                cc = corr(medfilt1(UpState1(end,:),(10/dt),'truncate')',medfilt1(UpState2(end,:),50,'truncate')');
                UPSTATECORRWINDOW = [UPSTATECORRWINDOW cc];
                
                %Calculating Firing Rate
                spikethr = -25;
                
                [~,spiketimesPY] = findpeaks(upstate1(j,preEventTime:preEventTime+updur1(j)),'MinPeakHeight',spikethr);
                FRPY = length(spiketimesPY)/((updur1(j))*(dt/1000));
                
                [~,spiketimesPV] = findpeaks(upstate2(j,preEventTime:preEventTime+updur1(j)),'MinPeakHeight',spikethr);
                FRPV = length(spiketimesPV)/((updur1(j))*(dt/1000));
                
                firingrate = [FRPY; FRPV];
                FIRINGRATE = [FIRINGRATE firingrate];
                
                %Firing Rate of PYstim, PYNon-Stim, PVStim, and PVNon-Stim
                if stimflag == 1
                    PreNonStimNumSpkPY = findpeaks(upstate1(j,preEventTime:stimtimes(1)),'MinPeakHeight',spikethr);
                    PreNonStimFRPY = length(PreNonStimNumSpkPY)/((stimtimes(1)-preEventTime)*(dt/1000));
                    
                    StimNumSpkPY = findpeaks(upstate1(j,stimtimes(1):stimtimes(2)),'MinPeakHeight',spikethr);
                    StimFRPY = length(StimNumSpkPY)/((stimtimes(2)-stimtimes(1))*(dt/1000));
                    
                    PostNonStimNumSpkPY = findpeaks(upstate1(j,stimtimes(2):preEventTime+3+updur1(j)),'MinPeakHeight',spikethr);
                    PostNonStimFRPY = length(PostNonStimNumSpkPY)/((updur1(j)-stimtimes(2)+3)*(dt/1000));
                    
                    PreNonStimNumSpkPV = findpeaks(upstate2(j,preEventTime:stimtimes(1)),'MinPeakHeight',spikethr);
                    PreNonStimFRPV = length(PreNonStimNumSpkPV)/((stimtimes(1)-preEventTime)*(dt/1000));                    
                    
                    StimNumSpkPV = findpeaks(upstate2(j,stimtimes(1):stimtimes(2)),'MinPeakHeight',spikethr);
                    StimFRPV = length(StimNumSpkPV)/((stimtimes(2)-stimtimes(1))*(dt/1000));
                    
                    PostNonStimNumSpkPV = findpeaks(upstate2(j,stimtimes(2):preEventTime+3+updur1(j)),'MinPeakHeight',spikethr);
                    PostNonStimFRPV = length(PostNonStimNumSpkPV)/((updur1(j)-stimtimes(2)+3)*(dt/1000));

                    stimfiringrate = [PreNonStimFRPY; StimFRPY; PostNonStimFRPY; PreNonStimFRPV; StimFRPV; PostNonStimFRPV];
                    STIMFIRINGRATE = [STIMFIRINGRATE stimfiringrate];
                    
                elseif stimflag == 0
                    PreVirtStimNumSpkPY = findpeaks(upstate1(j,(preEventTime:preEventTime+VStimLag)),'MinPeakHeight',spikethr);
                    PreVirtStimFRPY = length(PreVirtStimNumSpkPY)/((VStimLag)*(dt/1000)); 
                    
                    VirtStimNumSpkPY = findpeaks(upstate1(j,(preEventTime+VStimLag:preEventTime+VStimLag+VStimLen)),'MinPeakHeight',spikethr);
                    VirtStimFRPY = length(VirtStimNumSpkPY)/((VStimLen)*(dt/1000)); 
                    
                    PostVirtStimNumSpkPY = findpeaks(upstate1(j,(preEventTime+VStimLag+VStimLen:preEventTime+updur1(j))),'MinPeakHeight',spikethr);
                    PostirtStimFRPY = length(PostVirtStimNumSpkPY)/((VStimLag)*(dt/1000));
                    
                    PreVirtStimNumSpkPV = findpeaks(upstate2(j,(preEventTime:preEventTime+VStimLag)),'MinPeakHeight',spikethr);
                    PreVirtStimFRPV = length(PreVirtStimNumSpkPV)/((VStimLag)*(dt/1000)); 
                    
                    VirtStimNumSpkPV = findpeaks(upstate2(j,(preEventTime+VStimLag:preEventTime+VStimLag+VStimLen)),'MinPeakHeight',spikethr);
                    VirtStimFRPV = length(VirtStimNumSpkPV)/((VStimLen)*(dt/1000)); 
                    
                    PostVirtStimNumSpkPV = findpeaks(upstate2(j,(preEventTime+VStimLag+VStimLen:preEventTime+updur1(j))),'MinPeakHeight',spikethr);
                    PostirtStimFRPV = length(PostVirtStimNumSpkPV)/((VStimLag)*(dt/1000));                    
                    
                    viststimfiringrate = [PreVirtStimFRPY; VirtStimFRPY; PostirtStimFRPY; PreVirtStimFRPV; VirtStimFRPV; PostirtStimFRPV];
                    VIRTSTIMFIRINGRATE = [VIRTSTIMFIRINGRATE viststimfiringrate];
                end
                
                if graphics
                    figure(FiringRateFigure)
                    if isempty(FiringRateFigure.Children)
                        subplot(3,1,1)
                        plot((1:length(upstate1(j,:)))*dt,upstate1(j,:),'Color','k')
                        hold on
                        scatter((preEventTime+spiketimesPY)*dt,repmat(spikethr,1,length(spiketimesPY)),[],'b')
                        ylim([min(tracePY) 20])
                        xlim([0 length(upstate1)*dt])
                        subplot(3,1,2)
                        plot((1:length(upstate1(j,:)))*dt,upstate2(j,:),'Color','r')
                        hold on
                        scatter((preEventTime+spiketimesPV)*dt,repmat(spikethr,1,length(spiketimesPV)),[],'b')
                        ylim([min(tracePV) 20])
                        xlim([0 length(upstate1)*dt])
                        hold off
                        subplot(3,1,3)
                        plot((1:length(upstate1(j,:)))*dt,upstate3(j,:))
                        xlim([0 length(upstate1)*dt])
                        if stimflag == 1
                            hold on
                            scatter(stimtimes.*dt,[100 100],[],'k','MarkerFaceColor','k')
                            hold off
                        elseif stimflag == 0
                            hold on
                            scatter([],[],[],'k','MarkerFaceColor','k')
                            ylim([0 100])
                            hold off
                        end
                    elseif ~isempty(get(gca,'Children'))
                        %                     FiringRateFigure.Children(3).Children(2).XData = 0:updur1(j);
                        FiringRateFigure.Children(3).Children(2).YData = upstate1(j,:);
                        FiringRateFigure.Children(3).Children(1).XData = (preEventTime+spiketimesPY)*dt;
                        FiringRateFigure.Children(3).Children(1).YData = repmat(spikethr,1,length(spiketimesPY));
                        
                        %                     FiringRateFigure.Children(2).Children(2).XData = 0:updur1(j);
                        FiringRateFigure.Children(2).Children(2).YData = upstate2(j,:);
                        FiringRateFigure.Children(2).Children(1).XData = (preEventTime+spiketimesPV)*dt;
                        FiringRateFigure.Children(2).Children(1).YData = repmat(spikethr,1,length(spiketimesPV));
                        
                        
                        if stimflag == 1
                            FiringRateFigure.Children(1).Children(1).XData = stimtimes.*dt;
                            FiringRateFigure.Children(1).Children(2).YData = upstate3;
                        elseif stimflag == 0
                            FiringRateFigure.Children(1).Children(1).XData = [0 0];
                            FiringRateFigure.Children(1).Children(2).YData = upstate3;
                        end
                    end
                    drawnow
                end
                
                %Cross Correlation
                restPY = thrPY-thr;
                restPV = thrPV-thr;
                
                FiltUpStatePY = medfilt1(upstate1(j,:),round(5/dt),'truncate');
                FiltUpStatePV = medfilt1(upstate2(j,:),round(5/dt),'truncate');
                DataPY = FiltUpStatePY-mean(upstate1(j,:));
                DataPV = FiltUpStatePV-mean(upstate2(j,:));
                corrwindow50 = preEventTime-25/dt:preEventTime+25/dt;
                DataPY50 = DataPY(corrwindow50);
                DataPV50 = DataPV(corrwindow50);
                [C50,lag50] = xcorr(DataPY50,DataPV50);
                [C100,lag100] = xcorr(DataPY(250:750),DataPV(250:750));
                
                [~,peakindex50] = max(C50); peakshift50 = lag50(peakindex50)*dt;
                [~,peakindex100] = max(C100); peakshift100 = lag100(peakindex100)*dt;
                PeakShifts = [peakshift50; peakshift100];
                PEAKSHIFTS = [PEAKSHIFTS PeakShifts];
                
                center50 = find(lag50==0);
                center100 = find(lag100==0);
                
                AreaL50 = sum(C50(1:(center50-1)));
                AreaR50 = sum(C50((center50+1):end));
                area50 = [AreaL50; AreaR50];
                AREA50 = [AREA50 area50];
                
                AreaL100 = sum(C100(1:(center100-1)));
                AreaR100 = sum(C100((center100+1):end));
                area100 = [AreaL100; AreaR100];
                AREA100 = [AREA100 area100];
                
                if graphics
                    figure(XCorrFigure)
                    if isempty(XCorrFigure.Children)
                        subplot(3,1,1)
                        plot(lag50,C50,'LineWidth',2.0)
                        hold on
                        scatter(lag50(peakindex50),C50(peakindex50),[],'r')
                        line([0 0],[0 max(C50)],'Color',[0 0 0 0.3],'LineStyle','--')
                        title(['Pair# =' num2str(file) ', Cond=' num2str(cond) ', Trace#=' num2str(i) ', Upstate# in Trace=' num2str(j) ', Stim?' num2str(stimflag)])
                        hold off
                        subplot(3,1,2)
                        plot(lag100,C100,'LineWidth',2.0)
                        hold on
                        scatter(lag100(peakindex100),C100(peakindex100),[],'r')
                        line([0 0],[0 max(C100)],'Color',[0 0 0 0.3],'LineStyle','--')
                        hold off
                        subplot(3,1,3)
                        plot(DataPY50,'k','linewidth',2.0)
                        hold on
                        plot(DataPV50,'r','linewidth',2.0)
                    elseif ~isempty(get(gca,'Children'))
                        XCorrFigure.Children(3).Title.String = ['Pair# =' num2str(file) ', Cond=' num2str(cond) ', Trace#=' num2str(i) ', Upstate# in Trace=' num2str(j) ', Stim?' num2str(stimflag)];
                        XCorrFigure.Children(3).Children(3).XData = lag50;
                        XCorrFigure.Children(3).Children(3).YData = C50;
                        XCorrFigure.Children(3).Children(2).XData = lag50(peakindex50);
                        XCorrFigure.Children(3).Children(2).YData = C50(peakindex50);
                        XCorrFigure.Children(3).Children(1).YData = [0 max(C50)];
                        
                        XCorrFigure.Children(2).Children(3).XData = lag100;
                        XCorrFigure.Children(2).Children(3).YData = C100;
                        XCorrFigure.Children(2).Children(2).XData = lag100(peakindex100);
                        XCorrFigure.Children(2).Children(2).YData = C100(peakindex100);
                        XCorrFigure.Children(2).Children(1).YData = [0 max(C100)];
                        
                        XCorrFigure.Children(1).Children(2).YData = DataPY50;
                        XCorrFigure.Children(1).Children(1).YData = DataPV50;
                    end
                    drawnow
                    waitforbuttonpress
                end
            end
            
            clear upstart1 upstart2 startdiff upstate1 upstate2 upstate3 stimflag
        end
    end
    
    UPSTATES(1,:,:) = UpState1;
    UPSTATES(2,:,:) = UpState2;
    UPSTATES(3,:,:) = UpState3;
    
    STIMFLAG = logical(STIMFLAG);
    
    if ~exist('lag')
        lag = [];
    end
    
    META(file).StartDiff = STARTDIFF;
    META(file).MedianStartDiff = median(STARTDIFF);
    META(file).numEvents = length(STARTDIFF);
    META(file).UpStates = UPSTATES;
    META(file).UpStateCorr = UPSTATECORR;
    META(file).UpStateCorrWindow = UPSTATECORRWINDOW;
    META(file).MeanUpStateCorr = mean(UPSTATECORR);
    META(file).UpDurPY = UPDUR1;
    META(file).StimUpDurPY = UPDUR1(STIMFLAG);
    META(file).NonStimUpDurPY = UPDUR1(~STIMFLAG);
    META(file).StimFlag = STIMFLAG;
    META(file).Connectivity = Connectivity;
    META(file).FiringRate = FIRINGRATE;
    META(file).MeanFR = mean(FIRINGRATE,2);
    META(file).DistCells = str2double(DISTCELLS);
    META(file).PeakShifts = PEAKSHIFTS;
    META(file).MeanPeakShifts = mean(PEAKSHIFTS,2);
    META(file).Area50 = AREA50;
    META(file).MedianArea50 = median(AREA50,2);
    META(file).Area100 = AREA100;
    META(file).MedianArea100 = median(AREA100,2);
    META(file).TransType = TRANSTYPE;
    META(file).StimFiringRate = STIMFIRINGRATE; %IN HZ The order in rows is NonStimPY, StimPY, NonStimPV, StimPV
    META(file).MedianStimFiringRate = median(STIMFIRINGRATE,2); %IN HZ The order in rows is NonStimPY, StimPY, NonStimPV, StimPV
    META(file).MedianVirtStimFiringRate = median(VIRTSTIMFIRINGRATE,2);
    META(file).VirtStimFiringRate = VIRTSTIMFIRINGRATE;
    META(file).StimTimes = STIMTIMES;
end
toc
if 0
    for i=1:length(META)
        
        if META(i).numEvents>2
            clf
            subplot(3,1,1);
            imagesc(squeeze(META(i).UpStates(1,:,:)),[-80 -30])
            subplot(3,1,2);
            imagesc(squeeze(META(i).UpStates(2,:,:)),[-80 -30])
            subplot(3,1,3)
            meanPY = mean(squeeze(META(i).UpStates(1,:,:)));
            meanPV = mean(squeeze(META(i).UpStates(2,:,:)));
            plot(meanPY(1:2500))
            hold on
            plot(meanPV(1:2500))
            waitforbuttonpress
        end
        
        if META(i).numEvents>9
            DATA(i,:) = mean(META(i).UpStateCorr);
            DATA = DATA(DATA~=0)
        end
    end
    
    data=[META(ind).MeanStartDiff]
    [h p s c]=ttest(data)
    
    %% Setting Up Bullian Indexes
    ConMap = bin2dec(vertcat(META.Connectivity));
    
    eventind = vertcat(META.numEvents)>4; %More than 4 events
    cind = ConMap > 0; %Connected Cells
    noind = ConMap==0; %Non-Connected Cells
    pyrind = ConMap==2; %Pyr->PV
    pvind = or(ConMap==1,ConMap==3); %PV->Pyr
    recind = ConMap==3; %Reciprocal
    
    
    %Proper Indexes
    ceventind = and(cind,eventind);
    noeventind = and(noind,eventind);
    pyreventind = and(pyrind,eventind);
    pveventind = and(pvind,eventind);
    receventind = and(recind,eventind);
    
    %% Correlations of Upstates
    %Getting Proper Data and Transforming it
    CTransCorr = atanh([META(ceventind).MeanUpStateCorr]);
    NCTransCorr = atanh([META(noeventind).MeanUpStateCorr]);
    PYTransCorr = atanh([META(pyreventind).MeanUpStateCorr]);
    PVTransCorr = atanh([META(pveventind).MeanUpStateCorr]);
    ReciTransCorr = atanh([META(receventind).MeanUpStateCorr]);
    
    %Taking Averages
    MeanConn = mean(CTransCorr);
    MeanNotConn = mean(NCTransCorr);
    MeanPyrConn = mean(PYTransCorr);
    MeanPVConn = mean(PVTransCorr);
    MeanRecConn = mean(ReciTransCorr);
    
    %Standard Error
    SEConn = std(CTransCorr)/(sqrt(length(CTransCorr)));
    SENotConn = std(NCTransCorr)/(sqrt(length(NCTransCorr)));
    SEPyr = std(PYTransCorr)/(sqrt(length(PYTransCorr)));
    SEPV = std(PVTransCorr)/(sqrt(length(PVTransCorr)));
    SEREC = std(ReciTransCorr)/(sqrt(length(ReciTransCorr)));
    
    figure
    %     [hBar,hErrorbar] = barwitherr([SEConn SENotConn SEPyr SEPV],[MeanConn MeanNotConn MeanPyrConn MeanPVConn],'FaceColor','flat');
    [hBar,hErrorbar] = barwitherr([SENotConn SEREC],[MeanNotConn MeanRecConn],'FaceColor','flat');
    %     set(gca,'XTickLabel',{'Connected','Not Connected','PYR -> PV','PV -> PYR'})
    set(gca,'XTickLabel',{'Not Connected','Reciprocal'})
    set(hErrorbar,'Color','r','LineWidth',2)
    title('Transformed Correlation of Upstates')
    
    [h,p,s,c] = ttest2(ReciTransCorr,NCTransCorr)
    
    %% Relation of Distance to Connectivity
    
    DistConnect = [META(ceventind).DistCells];
    DistNonConn = [META(noeventind).DistCells];
    DistRecip = [META(receventind).DistCells];
    
    MeanDistConn = mean(DistConnect);
    MeanDistNonConn = mean(DistNonConn);
    MeanDistRec = mean(DistRecip);
    
    SEDCon = std(DistConnect)/(sqrt(length(DistConnect)));
    SEDNon = std(DistNonConn)/(sqrt(length(DistNonConn)));
    SEDRec = std(DistRecip)/(sqrt(length(DistRecip)));
    
    figure
    [hBar,hErrorbar] = barwitherr([SEDCon SEDNon SEDRec],[MeanDistConn MeanDistNonConn MeanDistRec],'FaceColor','flat');
    set(gca,'XTickLabel',{'Connected','Not-Connected','Reciprocal'})
    set(hErrorbar,'Color','r','LineWidth',2)
    title('Distance between the two cells')
    
    [h1,p1,s1,c1] = ttest2(DistRecip,DistNonConn)
    
    
    data = [DistConnect DistNonConn DistRecip];
    G    = [ones(1,length(DistConnect)) 2*ones(1,length(DistNonConn)) 3*ones(1,length(DistRecip))]
    boxplot(data,G);
    set(gca,'XTickLabel',{'Connected','Not-Connected','Reciprocal'})
    title('Distance between the two cells')
    
    
    
    %% Firing Rate
    %Getting Data
    AllFRData = [META(eventind).MeanFR];
    MeanFR = mean(AllFRData,2)';
    SEFR = (std(AllFRData,0,2)./sqrt(length(AllFRData)))';
    
    figure
    %     [hBar,hErrorbar] = barwitherr(SEFR,MeanFR,'FaceColor','flat');
    data = [AllFRData(1,:) AllFRData(2,:)];
    G    = [ones(1,length(AllFRData(1,:))) 2*ones(1,length(AllFRData(2,:)))]
    boxplot(data,G);
    set(gca,'XTickLabel',{'PYR','PV'},'FontSize',12)
    set(hErrorbar,'Color','r','LineWidth',2)
    title('Firing Rates','FontSize',15)
    ylabel('Firing Rate','FontSize',12)
    box off
    
    [h2,p2,s2,c2] = ttest2(AllFRData(1,:),AllFRData(2,:))
    
    %% Difference in Start of UpState
    StartData = ([META(eventind).MeanStartDiff]).*dt;
    MeanStartData = mean(StartData);
    
    [h3,p3,s3,c3] = ttest(StartData)
    
    %% Duration of Upstate in Perturbed and Nonpertubed
    
    count = 0; NonStimDur = []; StimDur = [];
    for i=1:length(META)
        stimflag = logical(META(i).StimFlag)
        if length(find(stimflag==0))>4 & length(find(stimflag==1))>4
            count = count+1;
            NonStimDur(count) = mean((META(i).UpDurPY(~stimflag))*dt);
            SEMNonStimDur(count) = std((META(i).UpDurPY(~stimflag))*dt)./sqrt(length((META(i).UpDurPY(~stimflag))))
            StimDur(count) = mean((META(i).UpDurPY(stimflag))*dt);
            SEMStimDur(count) = std((META(i).UpDurPY(stimflag))*dt)./sqrt(length((META(i).UpDurPY(stimflag))))
        end
    end
    
    data = [NonStimDur StimDur];
    G    = [ones(1,length(NonStimDur)) 2*ones(1,length(StimDur))]
    boxplot(data,G);
    set(gca,'XTickLabel',{'Non Stimulated','Stimulated'},'FontSize',12)
    set(hErrorbar,'Color','r','LineWidth',2)
    title('Average Duration of Upstate','FontSize',15)
    ylabel('Duration (ms)','FontSize',12)
    box off
    
    %     NonStimDur(eventind)
    %     StimDur(eventind)
end
