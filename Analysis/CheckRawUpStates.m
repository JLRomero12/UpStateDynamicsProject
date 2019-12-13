clear
close all

%% Choose Files and Load up the data
% [path_name, file_base, file_fmt] = data_info;
cd('..\ConnectivityExperiments\Pyr-SST-Recordings\')
load('SSTTdTomatoMETA')

%% Plotting Data

for i = 1:length(MasterMETA)
    
    StimFlag    = MasterMETA(i).StimFlag==1;
    NonStimFlag = MasterMETA(i).StimFlag==0;
    
    dataStim = squeeze(MasterMETA(i).UpStates(2,StimFlag,:));
    dataNon  = squeeze(MasterMETA(i).UpStates(2,NonStimFlag,:));
    
    if size(dataNon,1)>4 && size(dataNon,2)>4 %&& size(dataStim,1)>4 && size(dataStim,2)>4
        figure
        plot(dataStim','Color','b','LineWidth',2.0)
        hold on
        plot(dataNon','Color','r','LineWidth',2.0)
        rectangle('Position',[600 -80 10000 110],'FaceColor','none','EdgeColor','k','LineStyle','--','LineWidth',2.0);
        ylim([-80 30])
        hold off
        
        figure
        subplot(2,1,1)
        plot(dataStim(:,1:1000)','Color','b','LineWidth',2.0)
        subplot(2,1,2)
        plot(dataNon(:,1:1000)','Color','r','LineWidth',2.0)
        title(['Cell# ' num2str(i)])
        
        StimFigure   = gobjects(size(dataStim,1),1);
        NoStimFigure = gobjects(size(dataNon,1),1);
        
%         for i1 = 1:size(dataStim,1)
%             StimFigure(i1) = figure('Name',['Upstate ' num2str(i1)],'WindowState','maximized');
%             plot(dataStim(i1,:),'LineWidth',2.0);
%             hold on
%             rectangle('Position',[600 -80 10000 110],'FaceColor','none','EdgeColor','k','LineStyle','--','LineWidth',2.0);
%             ylim([-80 30])
%         end
%         StimTab = figs2tabs(StimFigure);
%         StimTab.Name = 'Stimulated Upstates';
        
        for i2 = 1:size(dataNon,1)
            NoStimFigure(i2) = figure('Name',['Upstate ' num2str(i2)],'WindowState','maximized');
            plot(dataNon(i2,:),'LineWidth',2.0);
            hold on
            rectangle('Position',[600 -80 10000 110],'FaceColor','none','EdgeColor','k','LineStyle','--','LineWidth',2.0);
            ylim([-80 30])
        end
        NoStimTab = figs2tabs(NoStimFigure);
        NoStimTab.Name = 'Not Stimulated Upstates';
        
        key = get(gcf, 'CurrentKey');
        while ~strcmp(key, 'return')
            waitforbuttonpress
            key = get(gcf, 'CurrentKey');
        end
        
        close all
    end
end