clear
close all
warning('off','signal:findpeaks:largeMinPeakHeight')

cd('C:\Users\BuonoLab\Desktop\Juan Luis\MainStudy\SOM-Halo\DoubleTrans')
Files = dir('S*.mat');
numFiles = length(Files);
[~,I] = sort([Files.datenum]);
Files = Files(I);
META(1:numFiles) = struct();
dt = 0.5;

Param.thr = 6;                %Threshold from the resting of the cell (median of data)
Param.SegLen = 7000/dt;       %This is the maximun length that we expect the upstate to be. (10000ms = 10 sec)
Param.MaxIntEvent = 250/dt;   %Max difference between the ThrCross and still be considered the same even (100 = 50ms)
Param.MaxIntUpstate = 250/dt; %Max difference between the ThrCross and still be considered the same even (100 = 50ms)
Param.MinAboveTime = 250/dt;   %Time that the event has to be above the threshold (Within the event) to count as Upstate (1000 = 500ms)
Param.MinOnOffInt = 250/dt;    %Time that the event has to be above the threshold (from beginning to end) to count as Upstate (2500 = 500ms)
Param.preEventTime = 100;     %Time given before the start of the upstate
Param.dt = dt;

StimFreqPoints = (1000/50)/dt - 1;
StimDelayPoints = 250/dt;

CheckUpstates = figure('Name','CheckUpstates');

for file = 16 %numFiles
    
    clearvars -EXCEPT META Files numFiles file FILE Param UpStateFigure DifferenceFigure CheckUpstates AllTrialsFigure FiringRateFigure XCorrFigure ...
        graphics VStimLag VStimLen StimFreqPoints StimDelayPoints
    
    load(Files(file).name)
%     load('S3C1_0211_.mat')
    
    UpState1    = [];
    UpState2    = [];
    UpState3    = [];
    STIMFLAG    = [];
    STIMTIMES   = [];
    UPDUR1      = [];
    FIRINGRATE  = [];
    UPSTATECORR = [];
    
    if ~isempty(find(strcmp(CondOrder,'NEW_BMI_LIGHT_FF'),1))
        UpStateIndex = find(strcmp(CondOrder,'NEW_BMI_LIGHT_FF'));
    end
    
    for cond = UpStateIndex
        
        ch1 = CH1{cond};
        ch2 = CH2{cond};
        ch3 = CH3{cond};
        dt = round(SampleRate{cond},1);
        
        Above =  (find(ch3>10)');
        AboveDiff = find(diff([-999 (Above)])>10);
        if length(AboveDiff)==1
            stimtype = 'NoStim';
        elseif length(AboveDiff)>1
            StimLen = ((Above(AboveDiff(2)-1)-Above(1))*(dt/6)+2);
            
            if StimLen>4900 && StimLen<5100
                stimtype = 'Long';
            elseif StimLen>490 && StimLen<510
                stimtype = 'Short';
            end
        end
        
        for i = 1:size(ch1,1)
            %Find the Upstate
            tracePY = ch1(i,:);
            tracePV = ch2(i,:);
            traceSTIM = ch3(i,:);
            [~, ~, ~, ~, ~, upstate1, updur1, upstart1, upend1]...
                =MAKESEGMENTS_STIM(tracePY,Param.thr,Param.MaxIntEvent,Param.MaxIntUpstate,Param.SegLen,Param.MinAboveTime,Param.MinOnOffInt,Param.preEventTime,0,traceSTIM);
            
            UpState1 = [UpState1; upstate1];
            
            
            for j = 1:length(upstart1)
                %% Selecting Other Channels
                upstate2(j,:) = ch2(i,upstart1(j)-Param.preEventTime:upstart1(j)+Param.SegLen);
                upstate3(j,:) = ch3(i,upstart1(j)-Param.preEventTime:upstart1(j)+Param.SegLen);
                
                %% Finding Stimulation
                StimTimes = find(upstate3(j,:)>10);
                if ~isempty(StimTimes)
                    firststim = StimTimes(1);
                    if firststim>Param.preEventTime+StimDelayPoints*0.9 && firststim < Param.preEventTime + StimDelayPoints*0.9 + StimDelayPoints
                        stimflag = 1;
                        lastpoint = find(diff([-9e9 StimTimes 9e9])>StimFreqPoints);
                        lastpoint = lastpoint(2);
                        stimtimes = [firststim; StimTimes(lastpoint-1)];
                    elseif firststim<Param.preEventTime+StimDelayPoints*0.9 || (firststim > Param.preEventTime + StimDelayPoints*0.9 + StimDelayPoints && firststim > Param.preEventTime + StimDelayPoints*0.9 + StimDelayPoints*2)
                        stimflag = -1;
                        stimtimes = [0; 0];
                    else
                        stimflag = 0;
                        stimtimes = [0; 0];
                    end
                elseif isempty(StimTimes)
                    stimflag = 0;
                    stimtimes = [0; 0];
                end
                
                
                %% Checking SelectedUpstates
                TraceTime = (1:length(tracePY))*dt/1000;
                UpStateTime = (1:length(upstate1(j,:)))*dt/1000;
                Light = (find(ch3(i,:)>10))*dt/1000;
                y = -10*ones(1,length(Light));
                
                figure(CheckUpstates)
                subplot(3,1,1)
                hold on
                plot(TraceTime,tracePY,'b','LineWidth',2)
                scatter(Light,y,[],'c','.')
                line([TraceTime(1) TraceTime(end)],[mode(tracePY)+Param.thr mode(tracePY)+Param.thr],'color','k','linestyle',':')
                line([upstart1(j)*dt/1000 upstart1(j)*dt/1000],[min(tracePY) max(tracePY)],'color','r','LineWidth',1.5,'LineStyle',':')
                line([upend1(j)*dt/1000 upend1(j)*dt/1000],[min(tracePY) max(tracePY)],'color','r','LineWidth',1.5,'LineStyle',':')
                ylim([-80 0])
                title(['Cell#:' num2str(file) ' Frame#:' num2str(i)])
                hold off
                subplot(3,1,2)
                hold on
                plot(UpStateTime,upstate1(j,:),'k','LineWidth',2)
                line([UpStateTime(1) UpStateTime(end)],[mode(tracePY)+Param.thr mode(tracePY)+Param.thr],'color','k','linestyle',':')
                line([Param.preEventTime*dt/1000 Param.preEventTime*dt/1000],[min(upstate1(j,:)) max(upstate1(j,:))],'color','k','linestyle',':','linewidth',1.5)
                line([(Param.preEventTime+updur1(j))*dt/1000 (Param.preEventTime+updur1(j))*dt/1000],[min(upstate1(j,:)) max(upstate1(j,:))],'color','k','linestyle',':','linewidth',1.5)
                title(['Duration:' num2str(updur1(j)*dt)])
                xlim([UpStateTime(1) UpStateTime(end)])
                hold off
                subplot(3,1,3)
                hold on
                plot(UpStateTime,upstate2(j,:),'r','LineWidth',2)
                line([Param.preEventTime*dt/1000 Param.preEventTime*dt/1000],[min(upstate1(j,:)) max(upstate1(j,:))],'color','r','linestyle',':','linewidth',1.5)
                line([(Param.preEventTime+updur1(j))*dt/1000 (Param.preEventTime+updur1(j))*dt/1000],[min(upstate1(j,:)) max(upstate1(j,:))],'color','r','linestyle',':','linewidth',1.5)
                plot(UpStateTime,upstate3(j,:),'k','LineWidth',1)
                title(['Stimmed: ' num2str(stimflag)])
                xlim([UpStateTime(1) UpStateTime(end)])
                hold off
                [x, y, button] = ginput(1);
                if button == 3
                    resp = input('new Stimflag (-1,0(noStim),1(stim) or 55 to change duration->');
                    if resp == 55
                        [x1, y1, ~] = ginput(1);
                        [x2, y2, ~]= ginput(1);
                        updur1(j) = (x2-x1)*(1000/dt);
                        updur1(j)
                        
                        resp = input('new Stimflag (-1,0(noStim),1(stim)->');
                    else
                        stimflag = resp;
                    end
                end
                %waitforbuttonpress
                clf
                STIMTIMES = [STIMTIMES stimtimes];
                STIMFLAG = [STIMFLAG stimflag];
                
                %% Firing Up-states
                spikethr = -25;
                
                [~,spiketimesPY] = findpeaks(upstate1(j,Param.preEventTime:Param.preEventTime+updur1(j)),'MinPeakHeight',spikethr);
                FRPY = length(spiketimesPY)/((updur1(j))*(dt/1000));
                
                [~,spiketimesPV] = findpeaks(upstate2(j,Param.preEventTime:Param.preEventTime+updur1(j)),'MinPeakHeight',spikethr);
                FRPV = length(spiketimesPV)/((updur1(j))*(dt/1000));
                
                firingrate = [FRPY; FRPV];
                FIRINGRATE = [FIRINGRATE firingrate];
                
                %% UpState Correlation
                FiltData1 = medfilt1(upstate1(j,Param.preEventTime:Param.preEventTime+updur1(j)),50,'truncate');
                FiltData2 = medfilt1(upstate2(j,Param.preEventTime:Param.preEventTime+updur1(j)),50,'truncate');
                c = corr(FiltData1',FiltData2');
                
                UPSTATECORR = [UPSTATECORR c];
            end
            
            UPDUR1   = [UPDUR1 updur1];
            
            if exist('upstate2')
                UpState2 = [UpState2; upstate2];
                UpState3 = [UpState3; upstate3];
            end
            
            clear upstate1 upstate2 upstate3 updur1
        end
    end
    UPSTATES(1,:,:) = UpState1;
    UPSTATES(2,:,:) = UpState2;
    UPSTATES(3,:,:) = UpState3;
    
    META(file).UpStates = UPSTATES;
    META(file).UpDurPY = UPDUR1;
    META(file).StimFlag = STIMFLAG;
    META(file).StimTimes = STIMTIMES;
    META(file).FiringRate = FIRINGRATE;
    META(file).MeanFR = mean(FIRINGRATE,2);
    META(file).Corr = UPSTATECORR;
    META(file).MeanCorr = mean(UPSTATECORR);
    
    StimIndex = META(file).StimFlag==1;
    NonStimIndex = META(file).StimFlag==0;
    StimUp = squeeze(META(file).UpStates(1,StimIndex,:));
    NonStimUp = squeeze(META(file).UpStates(1,NonStimIndex,:));
    
    subplot(1,2,1)
    imagesc(NonStimUp)
    caxis([-80 -40])
    subplot(1,2,2)
    imagesc(StimUp)
    caxis([-80 -40])
    
    NonStimUpDur = META(file).UpDurPY(1,NonStimIndex,:);
    StimUpDur = META(file).UpDurPY(1,StimIndex,:);
    g = [ones(1,length(NonStimUpDur)) 2*ones(1,length(StimUpDur))];
    
    figure
    boxplot([NonStimUpDur StimUpDur],g)
    
    NonStimUpDur = NonStimUpDur(NonStimUpDur>400);
    StimUpDur = StimUpDur(StimUpDur>400);
    g = [ones(1,length(NonStimUpDur)) 2*ones(1,length(StimUpDur))];
    figure
    boxplot([NonStimUpDur StimUpDur],g)
    
    META(file).PassNonStimUpDur = median(NonStimUpDur);
    META(file).PassStimUpDur = median(StimUpDur);
    if ~exist('stimtype','var')
        stimtype = 'Long'; %input('Is Stim Long or Short','s');
    end 
    META(file).StimType = stimtype;
end

filename = ['Cell-' num2str(file)];
save(filename,'META','Param')
%Notes: S1C1_1203_.mat looks like a really crappy cell.