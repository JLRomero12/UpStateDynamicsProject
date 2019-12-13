clear
close all
%% %%%%% Pyramidal to PV %%%%%%%%%%%
load OptoPyrPVConnectivity.mat
META = CONN;
clear CONN

load PyrPVConnectivity.mat
META = [META CONN];

ConMap = bin2dec(vertcat(META.Connectivity));

cind = ConMap > 0; %Connected Cells
noind = ConMap==0; %Non-Connected Cells
pyrind = ConMap==2; %Pyr->PV
pvind = ConMap==1;
recind = ConMap==3; %Reciprocal

ConPercent = [length(noind(noind~=0)) length(pyrind(pyrind~=0)) length(pvind(pvind~=0)) length(recind(recind~=0))];
figure
pie(ConPercent)
legend('Not Connected','PYR->PV','PV->PYR','Reciprocal','Location','eastoutside')
A = findobj(gca,'FontWeight','normal');
set(A,'FontWeight','bold','FontSize',12)
title('Connectivity Pyr and PV','FontSize',16,'HorizontalAlignment','left')


%% %%%%% Pyramidal to SST %%%%%%%%%%%
load OptoPyrSSTConnectivity.mat
META = CONN;
clear CONN

load PyrSSTConnectivity.mat
META = [META CONN];

ConMap = bin2dec(vertcat(META.Connectivity));

cind = ConMap > 0; %Connected Cells
noind = ConMap==0; %Non-Connected Cells
pyrind = ConMap==2; %Pyr->SST
pvind = ConMap==1; %SST->Pyr
recind = ConMap==3; %Reciprocal

ConPercent = [length(noind(noind~=0)) length(pyrind(pyrind~=0)) length(pvind(pvind~=0)) length(recind(recind~=0))];
figure
pie(ConPercent)
legend('Not Connected','PYR->SST','Reciprocal','Location','eastoutside')
A = findobj(gca,'FontWeight','normal');
set(A,'FontWeight','bold','FontSize',12)
title('Connectivity Pyr and SST','FontSize',16,'HorizontalAlignment','left')