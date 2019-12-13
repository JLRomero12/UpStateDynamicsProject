close all
clear

cd('C:\Users\BuonoLab\Desktop\Juan Luis\MainStudy\PV-Halo')
load MasterMETA.mat

stimtype = 'Short';
dt = Param.dt;
preStimTime = 275/dt;
stimdelay = 250/dt;


if strcmp(stimtype,'Long')
    postStimTime = 5100/dt;
    lenStim = 5000;
elseif strcmp(stimtype,'Short')
    postStimTime = 600/dt;
    lenStim = 500;
end
%% Init Figure Parameters
ImageSCFigure = figure('Name','ImageSC');

RawPlotFigure = figure('Name','RawPlots');
SP1 = subplot(2,1,1);
SP2 = subplot(2,1,2);

%% Finding Cells with the Proper Stimulation Type (Long or Short)
StimTypeCell = struct2cell(MasterMETA);
StimTypeCell = squeeze(StimTypeCell(end,:,:)); %Assumes that StimType will be the last field in the structuce
StimTypeInd = strcmp(StimTypeCell,stimtype);

%% Preload Variables and Structs
META = MasterMETA(StimTypeInd);
NewMETA(1,length(META)) = struct();

for i = 1:length(META)
    if length(META(i).UpDurPY)>9
        %% Loading up right data
        
        StimIndex       = META(i).StimFlag==1;
        NonStimIndex    = META(i).StimFlag==0;
        
        CH1     = squeeze(META(i).UpStates(1,:,:));
        CH3     = squeeze(META(i).UpStates(3,:,:));
        
        UpDurPY = META(i).UpDurPY;
        chunk   = zeros(size(CH3,1),preStimTime+postStimTime);
        
        %% Getting chunk centered around the first stimulus or 250 ms from onset (in non stimulated)
        for j = 1:size(CH3,1)
            if UpDurPY(j)>stimdelay
                StimStart = find(CH3(j,1:preStimTime*3)>10,1);
                
                if ~isempty(StimStart) && StimStart>500
                    chunk(j,:) = CH1(j,StimStart-preStimTime:StimStart+postStimTime-1);
                elseif isempty(StimStart)
                    chunk(j,:) = CH1(j,1+((Param.preEventTime+stimdelay)-preStimTime):(Param.preEventTime+stimdelay)+postStimTime);
                end
            end
        end
        
        stimchunk = chunk(StimIndex,:);
        stimchunk = stimchunk(any(stimchunk,2),:);
        
        nonstimchunk = chunk(NonStimIndex,:);
        nonstimchunk = nonstimchunk(any(nonstimchunk,2),:);
        
        AllChunk(i,:) = mean(nonstimchunk)-mean(stimchunk);
        
        %Saving for each pair
        NewMETA(i).StimType = stimtype;
        NewMETA(i).StimChunk = stimchunk;
        NewMETA(i).NonStimChunk = nonstimchunk;
        
        
        %% Plotting the Raw UpState Traces
        
        if ~isempty(stimchunk) && ~isempty(nonstimchunk)
            X = (1:length(stimchunk))*dt;
            
            figure(RawPlotFigure)
            subplot(SP1)
            plot(repmat(X',1,size(nonstimchunk,1)),nonstimchunk','LineWidth',1.5)
            hold on
            line([preStimTime*dt preStimTime*dt],[-80 0],'Color','k','LineStyle',':','LineWidth',2.0)
            line([(preStimTime*dt)+lenStim (preStimTime*dt)+lenStim],[-80 0],'Color','k','LineStyle',':','LineWidth',2.0)
            hold off
            ylim([-85 0])
            title('NonStimulated Upstates')
            
            subplot(SP2)
            plot(repmat(X',1,size(stimchunk,1)),stimchunk','LineWidth',1.5)
            hold on
            line([preStimTime*dt preStimTime*dt],[-80 0],'Color','k','LineWidth',2.0)
            line([(preStimTime*dt)+lenStim (preStimTime*dt)+lenStim],[-80 0],'Color','k','LineWidth',2.0)
            hold off
            ylim([-85 0])
            title('Stimulated UpStates')
            
            figure(ImageSCFigure)
            imagesc(AllChunk(i,:))
            title(['Cell #' num2str(i)])
            waitforbuttonpress
            clf
        end
    end
end

AllChunk = AllChunk(any(AllChunk,2),:);
figure(ImageSCFigure)
imagesc(AllChunk,[-10 10])
