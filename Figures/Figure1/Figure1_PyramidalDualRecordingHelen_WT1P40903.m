clear
close all
%first option
% load('WT1_P4_0903_.mat')
% UpState1 = CH1{2}(7,:);
% UpState2 = CH2{2}(7,:);

%second option 
% load('WT1_P1_1129_.mat')
% UpState1 = CH1{2}(7,:);
% UpState2 = CH2{2}(7,:);

%third option
% load('WT1_P2_0905_.mat')
% UpState1 = CH1{2}(8,:);
% UpState2 = CH2{2}(8,:);


%4TH option
% load('WT1_P2_0913_.mat')
% UpState1 = CH1{2}(?,:);
% UpState2 = CH2{2}(?,:);

%5TH OPTION
load('WT1_P3_0828_.mat')
UpState1 = CH1{2}(7,:);
UpState2 = CH2{2}(7,:);

dt = SampleRate(2);
x = (1:length(UpState1))*dt;

figure(1); set(gcf,'Position',[50 50 1200 400],'color','w')
plot(x,UpState1,'Color',[0 0.7 0],'LineWidth',2)
hold on
line([41000 41000],[-40 -20],'Color','k','LineWidth',2)
line([36000 41000],[-40 -40],'Color','k','LineWidth',2)
text(37500,-42,'5sec','FontSize',12,'FontWeight','bold')
h1 = text(42000,-26.5,'20mV','FontSize',12,'FontWeight','bold');
set(h1,'Rotation',270);

xlim([x(1) x(end)+8000])
ylim([-70 0])
[d1, X] = histcounts(UpState1);
plot(x(end)+4000+d1,conv(X, [0.5 0.5], 'valid'),'Color',[0 0.7 0],'LineWidth',2)
% axis tight
axis off

 print -dtiff WT1_P3_0828_PYR1



figure(2); set(gcf,'Position',[50 50 1200 400],'color','w')
plot(x,UpState2,'Color',[0 0.5 0],'LineWidth',2)
hold on
xlim([x(1) x(end)+8000])
ylim([-70 0])
[d1, X] = histcounts(UpState2);
plot(x(end)+4000+d1,conv(X, [0.5 0.5], 'valid'),'Color',[0 0.5 0],'LineWidth',2)
% axis tight
axis off

 print -dtiff WT1_P3_0828_PYR2
