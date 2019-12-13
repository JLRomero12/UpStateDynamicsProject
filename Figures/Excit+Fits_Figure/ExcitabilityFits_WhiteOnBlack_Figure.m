% RUN Figure1 first
% DATA PROBABLY COMES FROM PAIRED PV-PYR and SST-PYR data for connectivity and
% spontaneous (NOT OPTO)

%add zeros column
StepDur = .25; %sec

figure
set(gcf,'position',[80  340  880  260],'color','w')

for celltype = 1:3
   
   switch celltype
      case 1
         Spikes = [zeros(size(SpikesPYR,1),1) SpikesPYR];
         [numCells maxTraces] = size(SpikesPYR);
         color = [0 0.75 0];
         text = 'Pyr';
         ymax = 50;
      case 2
         Spikes = [zeros(size(SpikesPV,1),1) SpikesPV];
         [numCells maxTraces] = size(SpikesPV);
         color = 'r';
         text = 'PV';
         ymax = 100;
      case 3
         Spikes = [zeros(size(SpikesSST,1),1) SpikesSST];
         [numCells maxTraces] = size(SpikesSST);
         color = 'c';
         text = 'SST';
         ymax = 80;
   end
   INTEN = [0 0.05:0.05:0.05*(maxTraces)];
   Spks = nan(numCells,length(INTEN));
   
   count = 1; Accept = [];
   for i=1:numCells
      num = sum(Spikes(i,:));
      if num>1
         Accept = [Accept i];
         Spks(count,:) = Spikes(i,:);
         count = count+1;
      end
   end
   
   Spks    = Spks(1:count-1,:);
   Spks    = Spks/StepDur; %Hz
   SpkMean = nanmean(Spks,1);
   SpkSEM  =  nanstd(Spks)/sqrt(size(Spks,1));
   
   SP(celltype) = subplot(1,3,celltype);
   hold on
   labelfontsize = 16;
   xlabel('Intensity (nA)','fontweight','bold','fontsize',labelfontsize)
   if celltype == 1
      ylabel('Firing rate','fontweight','bold','fontsize',labelfontsize)
   end
   box off
   
   plot(INTEN,Spks','color',[0.2 0.2 0.2])
   hold on
   
   errorbar(INTEN,SpkMean,SpkSEM,'linewidth',4,'color','w')
   %set(gca,'ylim',[0 50])
   
   %% FIT
   F = @(param,x) param(1)*max(0,x-param(2)); %Threshold Linear activation function for Units ("Input-Output" function);
   X = INTEN;
   Y = SpkMean;
   
   param0 = [0, 0 ];
   param  = lsqcurvefit(F,param0,X,Y)
   Param(celltype,:) = param;
   
   %% Fit and Plot per cells
   clear ParamPerCell
   for i = 1:size(Spks,1)
      Y = Spks(i,:);
      Y = Y(Y>=0);
      param0 = [0, 0 ];
      ParamPerCell(i,:) = lsqcurvefit(F,param0,X(1:length(Y)),Y,[0 0],[inf inf]);
      yhat = F(ParamPerCell(i,:),X);
      plot(X,yhat,'linewidth',1,'color',[0.5 0.5 0.55])
      
      
   end
   PARAMPERCELL{celltype} = ParamPerCell;

  

   
   %% Plot Fit of Mean Excitability
   X = X(1):0.01:X(end);
   Yhat = F(param,X);
   plot(X,Yhat,'linewidth',4,'color',color)
   
   set(gca,'linewidth',2,'fontweight','bold','fontsize',14)
   set(gca,'xlim',[0 INTEN(end)],'ylim',[0 ymax])
   str = sprintf('\\fontsize{14}%s: \\fontsize{13}g=%3d, \\theta=%3.2f',text,round(Param(celltype,1)),Param(celltype,2));
   title(str,'color','w')
   set(gca,'color','k')

end

h=findall(0, 'Type', 'axes');
set(h(:),'xcolor','w','linewidth',2,'ycolor','w')
set(gca,'color','k')

set(gcf,'inverthardcopy','off')
set(gcf,'paperpositionmode','auto','color','k')
print -djpeg100 ExcitFits_WhiteOnBlack.jpg -r300



%%% Figure of Raw Excitabilities

colororder = [0.5 0.5 0.5; 1 1 1];
figure
set(gcf,'position',[80  340  880  260],'color','k')

SP(1) = subplot(1,3,1);
x1 = (0:3500)';
zWT = (Single(sampleWT).Data([1  8 ],1500:5000)');
xMat = repmat(x1,1,size(zWT,2));
zMatWT = zWT;
set(gca,'colororder',colororder)
hold on
plot(xMat,zMatWT,'linewidth',2)
ylim([-75 30])
set(gca,'color','k')
box off
axis off

SP(2) = subplot(1,3,2);
zPV = (Single(samplePV).Data([1 8 ],1500:5000)');
zMatPV = zPV;
set(gca,'colororder',colororder)
hold on
plot(xMat,zMatPV,'linewidth',2)
ylim([-75 30])
set(gca,'color','k')
box off
axis off

SP(3) = subplot(1,3,3);
zSST = (SSTMETA(sampleSST).Data{2}([1  8 ],6000:9500)');
zMatSST = zSST;
set(gca,'colororder',colororder)
hold on
plot(xMat,zMatSST,'linewidth',2)
ylim([-75 30])
set(gca,'color','k')
box off
axis off

set(gcf,'inverthardcopy','off')
set(gcf,'paperpositionmode','auto','color','k')
print -djpeg100 ExcitPyPVSST_Traces_WhiteOnBlack.jpg -r300

