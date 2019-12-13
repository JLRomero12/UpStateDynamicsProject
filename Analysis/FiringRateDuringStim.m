clear
close all

%% Choose Files and Load up the data
% [path_name, file_base, file_fmt] = data_info;
path_name = '..\..\PV_SST_Opto_Study_19\PairedUpstateData';
file_base = {'2Trans-SSTHalo-PairedMETA.mat'};

cd(path_name)
load(file_base{1})

%% Overal Variables
SpkAmp = 30;
restCriteria = -55;
DataChannel = 2; %Channel 1 is Pyr and Channel 2 is inhibitory neuron
CheckWork = 0;

%% For Non-Stimulated UpStates
dt        = Param.dt;
StartStim = (250/dt)+Param.preEventTime;
EndStim   = (225/dt)+StartStim;

%% Can Switch between Cheta and Halo

IsHalo = strfind(path_name,'Halo'); %determine if it is Halo stimulation
if isempty(IsHalo)
    StimUpThrs = 5000;
elseif ~isempty(IsHalo)
    StimUpThrs = 11000;
end

MasterMETA = PairedMETA;

%% Main Loop
for i = 1:length(MasterMETA)
    WantedUpStates = MasterMETA(i).StimFlag==1|MasterMETA(i).StimFlag==0;
    if isempty(IsHalo)
        SizeInd  = MasterMETA(i).UpDurPY>1000; %500 min length for upstates to be taken into consideration.
        data1    = squeeze(MasterMETA(i).UpStates(DataChannel,SizeInd&WantedUpStates,:));
        data3    = squeeze(MasterMETA(i).UpStates(3,SizeInd&WantedUpStates,:));
        updur    = MasterMETA(i).UpDurPY(SizeInd&WantedUpStates);
        StimFlag = MasterMETA(i).StimFlag(SizeInd&WantedUpStates);
    elseif ~isempty(IsHalo)
        SizeInd  = MasterMETA(i).UpDurPY>1000; %500 min length for upstates to be taken into consideration.
        data1    = squeeze(MasterMETA(i).UpStates(DataChannel,SizeInd&WantedUpStates,:));
        data3    = squeeze(MasterMETA(i).UpStates(3,SizeInd&WantedUpStates,:));
        updur    = MasterMETA(i).UpDurPY(SizeInd&WantedUpStates);
        StimFlag = MasterMETA(i).StimFlag(WantedUpStates&SizeInd);
    end
    
    if size(data1,2)>2 && size(data1(StimFlag==1,:),1)>4 && size(data1(StimFlag==0,:),1)>4
        CheckingFigure = gobjects(size(data1,1),1);
        FiringRate     = zeros(1,size(data1,1));
        
        for j = 1:size(data1,1)
            %% Fit everything to Channel 3's Stimulation or where stimulation would have happened
            resting = mode(data1(j,:));
            spkthr = resting+SpkAmp;
            
            if resting<restCriteria
                if StimFlag(j) == 1
                    stimInd = find(data3(j,1:StimUpThrs)>10);
                    stimInd = stimInd(1):stimInd(end);
                elseif StimFlag(j) == 0
                    stimInd = StartStim:EndStim;
                end
                
                if ~isempty(IsHalo)
                    stimInd = stimInd(1):updur(j);
                end
                
                %% Finding Firing Rate of the Stimulated region
                FakeSpkIndex  = find(data1(j,stimInd)>spkthr);
                RealSpkIndex  = find(diff([-99 FakeSpkIndex])>1);
                RealSpkIndex  = FakeSpkIndex(RealSpkIndex);
                spikeNum      = length(RealSpkIndex);
                FiringRate(j) = spikeNum/((length(stimInd)*Param.dt)/1000);
                
                %% Checking Work
                if (CheckWork)
                    CheckingFigure(j) = figure('Name',['Upstate ' num2str(j)],'WindowState','maximized');
                    plot(data1(j,:),'Color','k','LineWidth',2.0)
                    hold on
                    plot(data3(j,:),'Color','b','LineWidth',2.0)
                    rectangle('Position',[stimInd(1) -80 stimInd(end)-stimInd(1) 110])
                    scatter(RealSpkIndex+stimInd(1),repmat(spkthr,1,length(RealSpkIndex)),'filled','MarkerFaceColor','r')
                    line([1 length(data1(j,:))],[spkthr spkthr],'Color','k','LineStyle','--')
                    line([1 length(data1(j,:))],[resting resting],'Color','b','LineStyle','--')
                    hold off
                    title(['Number of Spikes =' num2str(spikeNum) ' Firing Rate =' num2str(FiringRate(j)) 'Cell # = ' num2str(i) 'StimFlag = ' num2str(StimFlag(j))])
                    ylim([-80 30])
                end
                
            else
                FiringRate(j) = NaN;
            end
        end
        %         CheckingTab = figs2tabs(CheckingFigure);
        %         CheckingTab.Name = 'CheckingFigure';
        %         addToolbarExplorationButtons(gcf)
        
        FR{1,i}       = FiringRate(StimFlag==0);  %Non-Stimulated On Top
        FR{2,i}       = FiringRate(StimFlag==1);  %Stimulated On the Bottom
        MedianFR(:,i) = [nanmedian(FiringRate(StimFlag==0));nanmedian(FiringRate(StimFlag==1))];
        MeanFR(:,i)   = [nanmean(FiringRate(StimFlag==0));nanmean(FiringRate(StimFlag==1))];
        
        if CheckWork
            key = get(gcf, 'CurrentKey');
            while ~strcmp(key, 'return')
                waitforbuttonpress
                key = get(gcf, 'CurrentKey');
            end
            close all
        end
    else
        FR{1,i}         = [];
        FR{2,i}         = [];
        MedianFR(:,i)   = [NaN;NaN];
        MeanFR(:,i)     = [NaN;NaN];
    end
end