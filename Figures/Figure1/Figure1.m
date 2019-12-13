clear
close all

ExamplePVCell = 'C2_0320_.mat';
UpStatePVTrace = [12 15 20];
GoodGreen = [0 0.7 0];

cd('..\..\ConnectivityExperiments\Pyr-PV-Recordings')
load(ExamplePVCell)

MainFigure = figure;
% MainFigure = figure('WindowState','Maximized');
set(gcf,'color','w');
ScreenSize = get(groot,'ScreenSize');
MainFigure.Position = [0 0 ScreenSize(3)/1.4 ScreenSize(4)/2]; %This is to Force the figure to have a vertical orientation

UpStateIndex = find(strcmp(CondOrder,'BMI_LIGHT'));
RawUpState1 = CH1(UpStateIndex);
RawUpState2 = CH2(UpStateIndex);
   
UpState1 = [];
UpState2 = [];

for i = 1:length(RawUpState1)
    upstate1 = RawUpState1{i};
    upstate2 = RawUpState2{i};
    UpState1 = [UpState1; upstate1];
    UpState2 = [UpState2; upstate2];    
end

dt = SampleRate(UpStateIndex(1));
x = (1:size(UpState1,2)*3)*dt;

figure(MainFigure)
plot(x,[UpState1(UpStatePVTrace(1),:) UpState1(UpStatePVTrace(2),:) UpState1(UpStatePVTrace(3),:)],'Color',GoodGreen,'LineWidth',1)
hold on
line([67250 67250],[-40 -10],'Color','k','LineWidth',2)
line([64250 67250],[-40 -40],'Color','k','LineWidth',2)
text(64450,-42.5,'3sec','FontSize',13,'FontWeight','bold','Color','k')
h1 = text(67998,-22,'30mV','FontSize',13,'FontWeight','bold','Color','k');
set(h1,'Rotation',270);
xlim([x(1) x(end)+500])
% axis tight
axis off

cd('..\..\Figures\Figure1')
print '-dtiffn' Figure1