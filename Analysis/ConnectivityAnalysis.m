clear
close all
warning('off','signal:findpeaks:largeMinPeakHeight')

cd('..\ConnectivityExperiments\Pyr-SST-Recordings')

Files = dir('S*.mat');
numFiles = length(Files);

% CONN(1:numFiles) = struct;

%% INIT GRAPHICS
ConnectivityFigure = figure;
CF1 = subplot('Position',[0.05 0.55 0.9 0.425]);
CF2 = subplot('Position',[0.05 0.05 0.9 0.425]);

%% Main Loop
for file = 14 %1:numFiles
    clearvars('ch1','ch2','PYRSig','PVSig','UpStateIndex')
    
    load(Files(file).name)
    CondOrder = CondOrder(~cellfun(@isempty,CondOrder));
    ConnectivityIndex = find(contains(CondOrder,'CONNECTIVITY'));
    EXCITIndex = find(strcmp(CondOrder,'EXCIT'));
    
    if ~isempty(Connectivity)
        
        CONN.Connectivity = Connectivity;
        CONN.Filename = Files(file).name;
        
    elseif isempty(Connectivity) && ~isempty(EXCITIndex)
        
        CondInten = [ExcitInten{:}];
        ExcitIndex = CondInten<99;
        CondInten = CondInten(ExcitIndex);
        
        ch1 = CH1(strcmp(CondOrder,'EXCIT'));
        ch1 = cell2mat(ch1');
        ch1 = ch1(ExcitIndex,:);
        ch2 = CH2(strcmp(CondOrder,'EXCIT'));
        ch2 = cell2mat(ch2');
        ch2 = ch2(ExcitIndex,:);
        
        for i = 3:size(ch1,1)
            
            figure(ConnectivityFigure)
            subplot(CF1)
            plot(ch1(i,:),'k','LineWidth',2.0)
            hold on
            [~,spLocPY] = findpeaks(ch2(i,6000:9000),'MinPeakHeight',max(ch2(i,6000:9000))-20);
            spLocPY = (spLocPY+6000);
            stimLinesPY = repmat(spLocPY,2,1);
            line(stimLinesPY,repmat([-85;0],1,length(stimLinesPY)),'Color',[1 0 0 0.5],'LineStyle','--','LineWidth',0.25)
            ylim([min(ch1(i,:))-5 max(ch1(i,:))+5])
            hold off
            
            subplot(CF2)
            plot(ch2(i,:),'r','LineWidth',2.0)
            hold on
            [~,spLocPV] = findpeaks(ch1(i,2000:4500),'MinPeakHeight',max(ch1(i,2000:4500))-20);
            spLocPV = (spLocPV+2000);
            stimLinesPV = repmat(spLocPV,2,1);
            line(stimLinesPV,repmat([-85;0],1,length(stimLinesPV)),'Color',[0 0 0 0.5],'LineStyle','--','LineWidth',0.25)
            ylim([min(ch2(i,:))-5 max(ch2(i,:))+5])
            hold off
            
            waitforbuttonpress
            
            %             key = get(gcf, 'CurrentKey');
            %             while ~strcmp(key, 'return')
            %                 waitforbuttonpress
            %                 key = get(gcf, 'CurrentKey');
            %             end
            
        end
        Connectivity = input('Connectivity?','s') ;
        CONN.Connectivity = Connectivity;
        CONN.Filename = Files(file).name;
        
    end
    
    filename = ['Cell-' num2str(file) '_' Files(file).name];
    save(filename,'CONN')
    
end


