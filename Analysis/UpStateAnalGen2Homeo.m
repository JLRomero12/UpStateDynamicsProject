clear
close all

cd('C:\Users\BuonoLab\Desktop\Juan Luis\HomeoPVCre\HomeoPVCompressedData')
Files = dir('P*.mat');
numFiles = length(Files);
META(1:numFiles) = struct();

CheckUpstates = figure('Name','CheckUpstates');

for file = 1 %numFiles
    
    clearvars -EXCEPT META Files numFiles file Param UpStateFigure DifferenceFigure CheckUpstates AllTrialsFigure FiringRateFigure XCorrFigure ...
        graphics VStimLag VStimLen StimFreqPoints StimDelayPoints
    
    load(Files(file).name)
    
    UpState1 = [];
    UPDUR1 = [];
    
    if ~isempty(find(strcmp(CondOrder,'SPONT'),1))
        UpStateIndex = find(strcmp(CondOrder,'SPONT'));
    end
    
    for cond = UpStateIndex
        
        dt = SampleRate{cond};
        Param.thr = 5;             %Threshold from the resting of the cell (median of data)
        Param.SegLen = 10000/dt;      %This is the maximun length that we expect the upstate to be.
        Param.MaxIntEvent = 25/dt;   %Max difference between the ThrCross and still be considered the same event
        Param.MaxIntUpstate = 100/dt; %Max difference between the ThrCross and still be considered the same upstate
        Param.MinAboveTime = 0/dt; %Time that the event has to be above the threshold (Within the event) to count as Upstate
        Param.MinOnOffInt = 100/dt;  %Time that the event has to be above the threshold (from beginning to end) to count as Upstate
        Param.preEventTime = 50;  %Time given before the start of the upstate
        StimFreqPoints = (1000/50)/dt - 1;
        StimDelayPoints = 250/dt;
        Param.dt = dt;
        
        ch1 = CH1{cond};
        dt = round(SampleRate{cond},1);
        
        for i = 1:size(ch1,1)
            %Find the Upstate
            tracePY = ch1(i,:);
            
            [~, ~, ~, ~, ~, upstate1, updur1, upstart1, upend1]...
                =MAKESEGMENTS(tracePY,Param.thr,Param.MaxIntEvent,Param.MaxIntUpstate,Param.SegLen,Param.MinAboveTime,Param.MinOnOffInt,Param.preEventTime,0);
            
            UpState1 = [UpState1; upstate1];
            
            
            for j = 1:length(upstart1)                                       
                %% Checking SelectedUpstates
                TraceTime = (1:length(tracePY))*dt/1000;
                UpStateTime = (1:length(upstate1(j,:)))*dt/1000;
                
                figure(CheckUpstates)
                hold on
                plot(TraceTime,tracePY,'b','LineWidth',2)
                line([TraceTime(1) TraceTime(end)],[mode(tracePY)+Param.thr mode(tracePY)+Param.thr],'color','k','linestyle',':')
                line([upstart1(j)*dt/1000 upstart1(j)*dt/1000],[min(tracePY) max(tracePY)],'color','r','LineWidth',1.5,'LineStyle',':')
                line([upend1(j)*dt/1000 upend1(j)*dt/1000],[min(tracePY) max(tracePY)],'color','r','LineWidth',1.5,'LineStyle',':')
                ylim([-80 0])
                title(['Cell#:' num2str(file) ' Frame#:' num2str(i)])
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
                
            end
            
            UPDUR1   = [UPDUR1 updur1];
                      
            clear upstate1 updur1
        end
    end
    UPSTATES = UpState1;
    
    META(file).UpStates = UPSTATES;
    META(file).UpDurPY = UPDUR1;
    
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
    
end

% save temp
%Notes: S1C1_1203_.mat looks like a really crappy cell.