clear
close all

load('MasterMETA')

MinDuration = 250;
TimLocPad = 600;
SuccLength = 750;

First3 = zeros(length(MasterMETA),3);
Last3  = zeros(length(MasterMETA),3);

for i = 1:length(MasterMETA)
    if length(MasterMETA(i).UpDurPY)>9
        StimIndex       = MasterMETA(i).StimFlag==1;
        NonStimIndex    = MasterMETA(i).StimFlag==0;
        
        NonStimUpDur    = (MasterMETA(i).UpDurPY(1,NonStimIndex,:))*Param.dt;
        StimUpDur       = (MasterMETA(i).UpDurPY(1,StimIndex,:))*Param.dt;
        
        PassNon         = NonStimUpDur(NonStimUpDur>MinDuration);
        PassStim        = StimUpDur(StimUpDur>MinDuration);
        
        if length(PassNon)>5
            First3(i,:) = PassNon(1:3);
            Last3(i,:)  = PassNon(end-2:end);
        else
            First3(i,:) = [NaN,NaN,NaN];
            Last3(i,:)  = [NaN,NaN,NaN];
        end
        
        SuccIndex       = find(PassStim<(SuccLength*1.1));
        FailIndex       = find(PassStim>(SuccLength*1.1));
        
        NonStimUpStates = squeeze(MasterMETA(i).UpStates(1,NonStimIndex,:));
        NonStimUpStates = NonStimUpStates(NonStimUpDur>MinDuration,:);
        StimUpStates    = squeeze(MasterMETA(i).UpStates(1,StimIndex,:));
        StimUpStates    = StimUpStates(StimUpDur>MinDuration,:);
        
        SuccStimUp      = StimUpStates(SuccIndex,:);
        FailStimUp      = StimUpStates(FailIndex,:);
        SuccDur         = PassStim(SuccIndex);
        FailDur         = PassStim(FailIndex);
        
        SuccRatio       = length(SuccDur)/length(PassStim);
        
%         figure
%         subplot(2,1,1)
% %         plot(SuccStimUp')
%         imagesc
%         title('Successfull Stimulation')
%         subplot(2,1,2)
%         plot(FailStimUp')
%         title('Failed Stimulation')
% %         waitforbuttonpress
        
        %% TimeLocked Plotting
        StimVectors     = squeeze(MasterMETA(i).UpStates(3,StimIndex,:));
        StimVectors     = StimVectors(StimUpDur>MinDuration,:);
        for j = 1:size(StimVectors,1)
            StimPoint  = find(StimVectors(j,:)>10,1);
            TimLocStimUp(j,:)    = StimUpStates(j,StimPoint-TimLocPad:StimPoint+TimLocPad);
        end
        figure
        subplot(2,1,1)
        hold on
        plot(TimLocStimUp')
        line([TimLocPad,TimLocPad],[min(min(TimLocStimUp)) max(max(TimLocStimUp))],'Color','k','LineStyle',':','LineWidth',2.0)
        title('Time Locked All Stimulated Upstates')
        hold off
        subplot(2,1,2)
        hold on
        plot(mean(TimLocStimUp,1))
        line([TimLocPad,TimLocPad],[min(min(TimLocStimUp)) max(max(TimLocStimUp))],'Color','k','LineStyle',':','LineWidth',2.0)
        title('Time Locked Average of All Stimulated Upstates')
        hold off
               
%                 g = [ones(1,length(PassNon)) 2*ones(1,length(PassStim))];
%                 x1 = ones(1,length(PassNon));
%                 x2 = 2*ones(1,length(PassStim));
%         
%                 figure
%                 subplot(3,1,1)
%                 boxplot([PassNon PassStim],g)
%                 hold on
%                 scatter(xC:\Users\BuonoLab\Desktop\Juan Luis\MainStudy\PV-Cheta1,PassNon,[],'r')
%                 scatter(x2,PassStim,[],'r')
%                 hold off
%                 subplot(3,1,2)
%                 imagesc(Param.dt:Param.dt:Param.SegLen*Param.dt,[],NonStimUpStates,[-80 -40])
%                 subplot(3,1,3)
%                 imagesc(Param.dt:Param.dt:Param.SegLen*Param.dt,[],StimUpStates,[-80 -40])
 
                
                
        NewMETA(i).PassNonStimUpDur = median(PassNon);
        NewMETA(i).PassStimUpDur = median(PassStim);
        NewMETA(i).SuccStimUp = SuccStimUp;
        NewMETA(i).FailStimUp = FailStimUp;
        NewMETA(i).StimUpStates = StimUpStates;
        NewMETA(i).NonStimUpStates = NonStimUpStates;
        NewMETA(i).SuccDur = median(SuccDur);
        NewMETA(i).FailDur = median(FailDur);
        NewMETA(i).SuccRatio = SuccRatio;
    end
    clear TimeLocStimUp
end


NonDur  = [NewMETA(:).PassNonStimUpDur];
StimDur = [NewMETA(:).PassStimUpDur];
SuccDur = [NewMETA(:).SuccDur];
FailDur = [NewMETA(:).FailDur];

figure
boxplot([NonDur' StimDur'])

[p,h,t,] = signrank(NonDur,StimDur);

DeltaUpDur = NonDur-StimDur;
DeltaFirstLast = median(First3,2)-median(Last3,2);
[~,sortedIndex] = sort(DeltaUpDur,'descend');




% %%
%
% for i = 1:length(MasterMETA)
%     StimIndex = MasterMETA(i).StimFlag==1;
%     NonStimIndex = MasterMETA(i).StimFlag==0;
%
%     figure
%
%
%
% end
