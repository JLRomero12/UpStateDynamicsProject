%%% Spont.m
% looks at a number of different spontaneous parameters
% ThrCross
% AreaAboveThr
% TimeAboveThr
% MeanCorrelation between segments (based on events)
% Sample is the Matrix of events

clear
global MOUSEBUTTON            %%% USED IN SpontPat_RESIZEAXES
close all


n = 1;

% filename{n} = 'ENTR033009_1A_5s'; n=n+1;
% filename{n} = 'ENTR033009_1B_5s'; n=n+1;
% filename{n} = 'ENTR032809_1A_5s'; n=n+1;
% filename{n} = 'ENTR032809_1B_5s'; n=n+1;
% filename{n} = 'ENTR040109_1A_5s'; n=n+1;
% filename{n} = 'ENTR040109_1B_5s';
%filename{n} = 'ENTR040909_1A_5s'; n=n+1;
filename{n} = 'ENTR040909_1B_5s'; n=n+1;

NUMFILES = n;

%%% IMPORTANT PARAMETERS %%%
VIEWTHRESHOLDING = 0;        %interactive
LOOKCORRELATIONS = 1;        %CALCULATE SIMILARITY BETWEEN TRACES
VIEWCORRS = 1;               %interactive...see traces and correlation values
timestep = 1;                % ASSUMES TIMES STEP = 1 ms (sample 1 kHz)
dt = 0.2;                    %timestep before resampling (in ms)
minDur=500;                  %min event duration for event to be used in calculating IEIs
SegLength = 3000;            %length of sample after onset of event (in ms)
preEventTime= -100;  %-10    %amount of baseline to include in Sample segements
preEventTimeCor = -25;       %amount of time before Event to be used in calculating correlation
postEventTimeCor = -preEventTimeCor;   %time after event used in calculating correlation
MaxInt = 100;               %Max allowable interval to still be 1 event; note that if MaxInt=50 and events = 1 49 98 it will only select event = 1;
MinOnOffInt= 500;          % make sure the segment onset and offset > then (elimnate single events)
MinUpTime = 50;            %USED TO REMOVE EPSP (i.e., events that are not UP States)
CorrThr = 0.5  ;            %Criteria for a match
MEDIANTHR = 1;              %do not use absolute Thr, but median + constant 
Thr = 5;                    %mv above median to be used for theresholding
COUNT = 0;
TotalSegLength=SegLength-preEventTime+1;
%%%% Fourier Analysis Variables %%%%
Fs=1000;                   % Sampling frequency
T = 1/Fs;                   % Sample time
L = 117500;                % Length of signal
t = (0:L-1)*T;              % Time vector
NFFT = 2^nextpow2(L);       % Next power of 2 from length of y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for n=1:NUMFILES
    clear data Sample OFF count SlpEventTimes EventTimes CH1 CH2 Type
    DATA =[]; Sample=[]; SegDur=[]; IEIs=[]; Spec=[];
    TotalNumEvents = 0; TotalAreaAboveThr = 0; TotalTimeAboveThr=0; TotalThrCross=0; 
    clear thr;
    file = filename{n};
    fprintf('%2d: %s\n',n,file);
    eval(['load ' file])
    clf;
    %plot(CH1{1}');
    %title('ACCEPT?  L = YES   R = NO');
    %[X, Y, BUTTON]=ginput(1);
    %if BUTTON==1
    for i=1:length(Cond)
        if strcmp(char(Cond(i)),'SPONT')
            SpontCond=i;
        end        
    end
    
    COUNT=COUNT+1;
    DATA = [CH1{SpontCond}];
    [numTraces numpoints] = size(DATA); ResampleDATA=spalloc(numTraces,numpoints/(timestep/dt),1);
    %%% SMOOTH AND DECIMATE %%%
    for i=1:numTraces
        trace = DATA(i,:);
        %SMOOTHTRACE WITH 4 ms windows (+/- 2ms)
        trace=SMOOTHTRACE(trace,20); %20 time points 
        %Resample at 1 ms step
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        trace=trace(1:1/dt:length(trace));      %%%Resample to 1ms           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ResampleDATA(i,:)=trace;
    end   
    DATA = ResampleDATA;
    [numTraces numpoints] = size(DATA);   
    MaxTime=numpoints;          %%%%numpoints in each sweep in ms %%%            
    TotalTime = (numpoints*numTraces)/1000;    %%%%total time of all sweeps in s %%%%
    fprintf('numTraces=%4d numpoint=%6d\n',numTraces,numpoints);
    time = timestep:timestep:MaxTime;
    clf; hold off

    for i=1:numTraces
        cla; hold on; set(gca,'units','pixels');
        MOUSEBUTTON=1; %necessary to reset the tag changed in RESIZEAXES
        trace = DATA(i,:); 
        plot(time,trace)
        thr=median(trace)+Thr;
        line([0 MaxTime], [thr thr],'color','g')
        
        %%%%%%%%%%%% FIND EVENTS %%%%%%%%%%%%%%%
        [eventtimes, offtimes, tabove, areaabove, numThrX, sample, segdur]= ...
        MAKESEGMENTS(trace,Thr,MaxInt,SegLength,MinUpTime,MinOnOffInt,preEventTime);
        %%% PLOT AND FILL EVENT  %%%%%%%%%%%%%%%
        numevents=length(eventtimes);
        for j=1:numevents
            line([(eventtimes(j)) (eventtimes(j))],[-30 -10],'color',[1 0 0],'linewidth',[3]);
            segment=trace(eventtimes(j):offtimes(j));
            segment(segment<thr)=thr;
            fill([eventtimes(j):offtimes(j) offtimes(j)],[segment thr],[0.5 0.5 0.5]);
%             if j<numevents
%                 IEIs=[IEIs eventtimes(j+1:numevents)-eventtimes(j)];
%             end
        end
        
        ieis=eventtimes(offtimes-eventtimes>minDur);
        IEIs=[IEIs diff(ieis)];
        Sample=[Sample; sample];
        SegDur = [SegDur segdur];
        TotalNumEvents = TotalNumEvents + numevents;
        TotalTimeAboveThr=TotalTimeAboveThr + tabove;
        TotalAreaAboveThr=TotalAreaAboveThr + areaabove;
        TotalThrCross = TotalThrCross + numThrX;
        str = sprintf('Sweep #%d',i);
        texthand = text(100,0,str);
        set(gca,'ylim',[-90 20]);
        set(gca,'xlim',[time(1) time(length(time))]);
        drawnow
        %%% WAIT FOR RESPONSE BEFORE NEXT TRACE %%%
        if (VIEWTHRESHOLDING) 
            SpontPat_RESIZEAXES;  %THIS IS CALLED HERE TO RESET THE INITALAXES
            while MOUSEBUTTON ~=3
            drawnow;
            end
        end 
    end    %%%%%%end numTraces %%%%%
    
      
   
   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% TRY TO CLUSTER SAMPLE TRACES, FIGURE OUT WHICH OF THE MATCHES ARE PART A GROUP %%%
        %%% USED THE TRACES FROM STORETRACE (the one's you RIGHT clicked on)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
        MeanCorr=NaN; NumMatches=NaN; TotalPairs=NaN; SimilMatrix=NaN; MeanUpStates=NaN; NumUpStates=NaN;
        if (LOOKCORRELATIONS)
            [numSamples numSamplePoints]=size(Sample);
            [row col]=size(Sample);
            %%% MAKESURE AT LEAST THREE SAMPLES THAT SATISFIED minUpTime were found
            if (row>3)      
                clear SimilMatrix
                fprintf('CALCULATING CORRELATIONS BETWEEN TRACES\n')
                for i=1:numSamples-1
                    for j=i:numSamples
                        clf; hold on;
                        maxLen=max(SegDur(i),SegDur(j))+postEventTimeCor-preEventTime;
                        if (maxLen>numSamplePoints) 
                            maxLen=numSamplePoints; 
                        end
                        dummy=(-preEventTime+preEventTimeCor);
                        trace1=Sample(i,(dummy+1):(maxLen) ); 
                        trace2=Sample(j,(dummy+1):(maxLen) ); 
                        corr = corrcoef(trace1,trace2); corr=corr(1,2);
                        if (VIEWCORRS)
                            plot(trace1,'b'); plot(trace2,'r'); set(gca,'xlim',[0 TotalSegLength])
                            dummy=sprintf('corr=%6.3f',corr);
                            title(dummy)
                            drawnow
                            %waitforbuttonpress
                        end
                        SimilMatrix(i,j)=corr; SimilMatrix(j,i)=corr;
                    end
                end
                
         
                %%% CLUSTER %%%
                [row col]=size(SimilMatrix);
                if (row>3)
                    Y=pdist(SimilMatrix,'seuclid');
                    Z=linkage(Y,'average'); %average
                    ClusterVector=cluster(Z,0.01);
                    [clust clusterind] =sort(ClusterVector);
                    %%% ZERO DIAGONAL %%%
                    diagind=[0:numSamples-1]*numSamples+[1:numSamples];
                    SimilMatrix(diagind)=0;
         
                    %%% HOW MANY CORRELATIONS ABOVE AN ARBITRARY ThrESHOLD
                    CorrData=tril(SimilMatrix);
                    [dummy1 dummy2 CorrVector]=find(tril(SimilMatrix,-1));
                    MeanCorr=mean(CorrVector);
                    NumMatches=length(find(CorrData>CorrThr));
                    TotalPairs=numSamples*(numSamples-1)/2;
                end
         
                MeanUpStates=mean(SegDur);
                NumUpStates=length(SegDur);
            end
        end 
       
 
        %%% end LOOKCORRELATIONS %%%
   
        TotalAreaAboveThr=TotalAreaAboveThr/(1000*TotalTime);
        TotalNumEvents=TotalNumEvents/TotalTime;
        TotalTimeAboveThr=TotalTimeAboveThr/(1000*TotalTime);
        TotalThrCross=TotalThrCross/TotalTime;
   
        fprintf('MinUpTime=%4d, MinIEI=%4d, preEventTime=%4d, Thr=%5.2f, CorrThr=%4.2f\n',MinUpTime,MaxInt,preEventTime,Thr,CorrThr);
        fprintf('NumMat=%4d (%4d) : %5.3f | MeanCorr=%5.2f\n',NumMatches,TotalPairs,NumMatches/TotalPairs,MeanCorr);
        fprintf('MeanUpStates=%5.3fms (%4d) \n',MeanUpStates,NumUpStates);
        fprintf('TimeAboveThr=%5.3f/s | TotalNumEvents=%6.3f/s | TotalArea=%8.4fmv/s (TotalTime=%5d)\n',TotalTimeAboveThr,TotalNumEvents,TotalAreaAboveThr,TotalTime);
        fprintf('\n')
        META_AreaAboveThr(COUNT)=TotalAreaAboveThr;
        META_Div(COUNT)=DIV;
        META_Filename(COUNT,:)=filename{n};
        META_IEI{COUNT}=IEIs;
        META_MatchIndex(COUNT)=NumMatches/TotalPairs;
        META_MeanCorr(COUNT)=MeanCorr;
        META_MeanUpStates(COUNT)=MeanUpStates;
        META_MatchProb(COUNT)=NumMatches/TotalPairs;
        META_NumEvents(COUNT)=TotalNumEvents;
        META_NumUpStates(COUNT)=NumUpStates;
        META_NumSeg(COUNT)=length(SegDur);
        META_OFF{COUNT}=SegDur;
        META_UpStates(COUNT)=MeanCorr;
        META_Sample{COUNT}=Sample;
        META_SegDur{COUNT}=SegDur;
        META_SimilMatrix{COUNT}=SimilMatrix;
        META_ClusterIndex{COUNT}=clusterind;
        META_TimeAboveThr(COUNT)=TotalTimeAboveThr;
        META_ThrCross(COUNT)=TotalThrCross;
        META_TotalPairs(COUNT)=TotalPairs;
        META_TotalTime(COUNT)=TotalTime;
        META_Type(COUNT,:)=Entrained;
        
    %end
end

%-----------------  END NUM FILES -------------------------------%

%%% SUMMARY DATA %%%

%-----------------------------------------------------------------------------%      

META_MeanDur=META_TimeAboveThr./META_NumEvents;
      
T=find(META_Type=='Y');
C=find(META_Type=='N');

NT=META_NumEvents(T);
NC=META_NumEvents(C);

Mean=[mean(NT) mean(NC)];
Mean1=Mean;

SEM(1)=std((NT)/sqrt(length(NT)));
SEM(2)=std((NC)/sqrt(length(NC)));
Sem=[SEM(1) SEM(2)];
Sem1=Sem;

AT=META_AreaAboveThr(T);
AC=META_AreaAboveThr(C);

Mean=[mean(AT) mean(AC)];
Mean2=Mean;

SEM(1)=std((AT)/sqrt(length(AT)));
SEM(2)=std((AC)/sqrt(length(AC)));
Sem=[SEM(1) SEM(2)];
Sem2=Sem;

TT=META_TimeAboveThr(T);
TC=META_TimeAboveThr(C);

Mean7=[mean(TT) mean(TC)];

SEM(1)=std((TT)/sqrt(length(TT)));
SEM(2)=std((TC)/sqrt(length(TC)));
Sem7=[SEM(1) SEM(2)];


NUT=META_NumUpStates(T);
NUC=META_NumUpStates(C);

NUT(find(isnan(NUT)))=0;
NUC(find(isnan(NUC)))=0;

Mean=[mean(NUT) mean(NUC)];
Mean8=Mean;

SEM(1)=std((NUT)/sqrt(length(NUT)));
SEM(2)=std((NUC)/sqrt(length(NUC)));
Sem=[SEM(1) SEM(2)];
Sem8=Sem;

MUT=META_MeanUpStates(T);
MUC=META_MeanUpStates(C);

Mean=[nanmean(MUT) nanmean(MUC)];
Mean9=Mean;

SEM(1)=std((MUT(find(~isnan(MUT))))/sqrt(length(MUT(find(~isnan(MUT))))));
SEM(2)=std((MUC(find(~isnan(MUC))))/sqrt(length(MUC(find(~isnan(MUC))))));
Sem=[SEM(1) SEM(2)];
Sem9=Sem;

X=[1 2];


clf
set(gcf,'color','w'); 
set(gca,'linewidth',[3],'fontweight','bold','fontsize',[26]);
subplot 221
hold on
X1=[0.8 1.8];
X2=[1.2 2.2];

Y1=[Mean1(1) Mean2(1)];
Y2=[Mean1(2) Mean2(2)];
S1=[Sem1(1) Sem2(1)];
S2=[Sem1(2) Sem2(2)];

h1=bar(X1,Y1,0.4);
set(h1,'facecolor',[0.3 0.3 0.3])
h2=bar(X2,Y2,0.4);
set(h2,'facecolor',[0.7 0.7 0.7])
xlim([0.5 2.5])
set(gca,'xtick',[1:1:2])
set(gca,'xticklabel',['Freq';'Area'])
legend(['Stimulated';'Control   ']);
%set(gca,'position',[0.18 0.11 0.775 0.815])

for i=1:2
    line([X1(i) X1(i)],[Y1(i)-S1(i) Y1(i)+S1(i)],'color','k','linewidth',[3]);
    line([X2(i) X2(i)],[Y2(i)-S2(i) Y2(i)+S2(i)],'color','k','linewidth',[3]);
end

subplot 222
hold on
bar(X,Mean7)
xlim([0.5 2.5])
set(gca,'xtick',[1:1:2])
set(gca,'xticklabel',['Stim';'Cont'])
title('Time Above Thr')

line([1 1],[Mean7(1)-Sem7(1) Mean7(1)+Sem7(1)],'color','k','linewidth',3);
line([2 2],[Mean7(2)-Sem7(2) Mean7(2)+Sem7(2)],'color','k','linewidth',3);

subplot 223
hold on
bar(X,Mean8)
xlim([0.5 2.5])
set(gca,'xtick',[1:1:2])
set(gca,'xticklabel',['Stim';'Cont'])
title('Num Up States')

line([1 1],[Mean8(1)-Sem8(1) Mean8(1)+Sem8(1)],'color','k','linewidth',3);
line([2 2],[Mean8(2)-Sem8(2) Mean8(2)+Sem8(2)],'color','k','linewidth',3);


subplot 224
hold on
bar(X,Mean9)
xlim([0.5 2.5])
set(gca,'xtick',1:1:2)
set(gca,'xticklabel',['Stim';'Cont'])
title('Mean Up State Dur')

line([1 1],[Mean9(1)-Sem9(1) Mean9(1)+Sem9(1)],'color','k','linewidth',3);
line([2 2],[Mean9(2)-Sem9(2) Mean9(2)+Sem9(2)],'color','k','linewidth',3);


figure
set(gcf,'color','w')

M1=[mean(META_MeanCorr(T)) mean(META_MeanCorr(C))];
Sem1=std(META_MeanCorr(T))/sqrt(length(T));
Sem2=std(META_MeanCorr(C))/sqrt(length(C));
S1=[Sem1 Sem2];

subplot 221
hold on
bar(X,M1)
xlim([0.5 2.5])
set(gca,'xtick',1:1:2)
set(gca,'xticklabel',['Stim';'Cont'])
title('Mean Correlation')

line([1 1],[M1(1)-S1(1) M1(1)+S1(1)],'color','k','linewidth',3);
line([2 2],[M1(2)-S1(2) M1(2)+S1(2)],'color','k','linewidth',3);

M2=[mean(META_MatchIndex(T)) mean(META_MatchIndex(C))];
Sem1=std(META_MatchIndex(T))/sqrt(length(T));
Sem2=std(META_MatchIndex(C))/sqrt(length(C));
S2=[Sem1 Sem2];

subplot 222
hold on
bar(X,M2)
xlim([0.5 2.5])
set(gca,'xtick',1:1:2)
set(gca,'xticklabel',['Stim';'Cont'])
title('Match Index')

line([1 1],[M2(1)-S2(1) M2(1)+S2(1)],'color','k','linewidth',3);
line([2 2],[M2(2)-S2(2) M2(2)+S2(2)],'color','k','linewidth',3);

TSamples=[];
for i=1:length(T)
    TSamples=[TSamples; META_Sample{T(i)}];
end

CSamples=[];
for i=1:length(C)
    CSamples=[CSamples; META_Sample{C(i)}];
end

subplot('position',[0.04 0.05 0.45 0.3])
imagesc(TSamples,[-60 -30])
subplot('position',[0.04 0.35 0.45 0.15])
plot(TSamples')
hold on
plot(mean(TSamples),'linewidth',3,'color','k')
xlim([0 length(TSamples)])
ylim([min(min(TSamples)) max(max(TSamples))])
axis off
title('Trained')

subplot('position',[0.53 0.05 0.45 0.3])
imagesc(CSamples,[-60 -30])
subplot('position',[0.53 0.35 0.45 0.15])
plot(CSamples')
hold on
plot(mean(CSamples),'linewidth',3,'color','k')
xlim([0 length(CSamples)])
ylim([min(min(CSamples)) max(max(CSamples))])
axis off
title('Control')


TIEIs=[];
for i=1:length(T)
    TIEIs=[TIEIs META_IEI{T(i)}];
end

CIEIs=[];
for i=1:length(C)
    CIEIs=[CIEIs META_IEI{C(i)}];
end

x=50:100:20000;
TDist=hist(TIEIs,x);
CDist=hist(CIEIs,x);

figure
set(gcf,'color','w')
subplot 221
bar(x,TDist)
ylim([0 10])
title('IEIs:  Trained')

subplot 222
bar(x,CDist)
ylim([0 10])
title('IEIs:  Control')

subplot 223
hold on
plot(TDist/sum(TDist),'color',[0 0.7 0])
plot(CDist/sum(CDist),'color',[0.7 0 0.7])
set(gca,'xticklabel',0:5:20)

mTIEI=[];
for i=1:length(T)
    mTIEI=[mTIEI mean(META_IEI{T(i)})];
end

mCIEI=[];
for i=1:length(C)
    mCIEI=[mCIEI mean(META_IEI{C(i)})];
end

M1=[mean(mTIEI) mean(mCIEI)];
Sem1=std(mTIEI)/sqrt(length(mTIEI));
Sem2=std(mCIEI)/sqrt(length(mCIEI));
S1=[Sem1 Sem2];

subplot 224
bar(X,M1)
xlim([0.5 2.5])
set(gca,'xtick',1:1:2)
set(gca,'xticklabel',['Stim';'Cont'])
title('Mean IEI')

line([1 1],[M1(1)-S1(1) M1(1)+S1(1)],'color','k','linewidth',3);
line([2 2],[M1(2)-S1(2) M1(2)+S1(2)],'color','k','linewidth',3);


%%%

TUpDurs=[];
for i=1:length(T)
    TUpDurs=[TUpDurs META_SegDur{T(i)}];
end

CUpDurs=[];
for i=1:length(C)
    CUpDurs=[CUpDurs META_SegDur{C(i)}];
end

x=50:100:20000;
TDurDist=hist(TUpDurs,x);
CDurDist=hist(CUpDurs,x);

figure
set(gcf,'color','w')
hold on
plot(x,TDurDist/sum(TDurDist),'color',[0 0.7 0],'linewidth',3)
plot(x,CDurDist/sum(CDurDist),'color',[0.7 0 0.7],'linewidth',3)
xlim([0 5000])
xlabel('Up State Duration')
ylabel('% Up States')

[X2, prob, df]=CHI2TEST([TDurDist(6:30); CDurDist(6:30)])