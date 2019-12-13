close all
clear

cd('C:\Users\BuonoLab\Desktop\Juan Luis\MainStudy\PV-Cheta')
load MasterMETA.mat

dt = Param.dt;
chunksize = 275/dt;
minlength = 250/dt;
stimdelay = 250/dt;
NewMETA(1,length(MasterMETA)) = struct();

%% Init Figure Parameters
ImageSCFigure = figure('Name','ImageSC');

RawPlotFigure = figure('Name','RawPlots');
SP1 = subplot(2,1,1);
SP2 = subplot(2,1,2);

for i = 1:length(MasterMETA)   %%ALERT!! CELLS 11 AND 12 ARE CHETA POSITIVE SO THEY ARE BAD FOR THIS ANALYSIS
    if length(MasterMETA(i).UpDurPY)>9 %&& i~=10 && i~=11
        %% Loading up right data 
        
        StimIndex       = MasterMETA(i).StimFlag==1;
        NonStimIndex    = MasterMETA(i).StimFlag==0;
        
        CH1 = squeeze(MasterMETA(i).UpStates(1,:,:));
        CH3 = squeeze(MasterMETA(i).UpStates(3,:,:));
        
        UpDurPY = MasterMETA(i).UpDurPY;
        chunk = zeros(size(CH3,1),chunksize*3.5);
        
        %% Getting chunk centered around the first stimulus or 250 ms from onset (in non stimulated)
        for j = 1:size(CH3,1)
            if UpDurPY(j)>minlength
                Ind = find(CH3(j,1:chunksize*3)>10,1);
                if ~isempty(Ind) && Ind>500
                    chunk(j,:) = CH1(j,Ind-chunksize:Ind+(chunksize*2.5)-1);
                elseif isempty(Ind)
                    chunk(j,:) = CH1(j,1+((Param.preEventTime+stimdelay)-chunksize):(Param.preEventTime+stimdelay)+(chunksize*2.5));
                end
            end
        end
        
        stimchunk = chunk(StimIndex,:);
        stimchunk = stimchunk(any(stimchunk,2),:);
        
        nonstimchunk = chunk(NonStimIndex,:);
        nonstimchunk = nonstimchunk(any(nonstimchunk,2),:);
        
        AllChunk(i,:) = mean(nonstimchunk)-mean(stimchunk);
        
        %Saving for each pair
        NewMETA(i).StimChunk = stimchunk;
        NewMETA(i).NonStimChunk = nonstimchunk;
        
        
        %% Plotting the Raw UpState Traces
        X = (1:length(stimchunk))*dt;
        
        figure(RawPlotFigure)
        subplot(SP1)
        plot(repmat(X',1,size(nonstimchunk,1)),nonstimchunk','LineWidth',1.5)
        hold on
        line([chunksize*dt chunksize*dt],[-80 0],'Color','k','LineStyle',':','LineWidth',2.0)
        line([(chunksize*dt)+500 (chunksize*dt)+500],[-80 0],'Color','k','LineStyle',':','LineWidth',2.0)
        hold off
        ylim([-80 0])
        title('NonStimulated Upstates')
        
        subplot(SP2)
        plot(repmat(X',1,size(stimchunk,1)),stimchunk','LineWidth',1.5)
        hold on
        line([chunksize*dt chunksize*dt],[-80 0],'Color','k','LineWidth',2.0)
        line([(chunksize*dt)+500 (chunksize*dt)+500],[-80 0],'Color','k','LineWidth',2.0)
        hold off
        ylim([-80 0])
        title('Stimulated UpStates')
        
        figure(ImageSCFigure)
        imagesc(AllChunk(i,:))
        title(['Cell #' num2str(i)])
%         waitforbuttonpress
        clf
    end
end

AllChunk = AllChunk(any(AllChunk,2),:);
figure(ImageSCFigure)
imagesc(AllChunk,[-10 10])
