function [eventtimes, offtimes, Tabove, areaabove, numThrX, upstate, updur, upstart, upend] ...
    =MAKESEGMENTS_STIM(trace,thr,MaxIntEvent,MaxIntUpstate,SegLen,MinAboveTime,MinOnOffInt,preEventTime,KurtosisFlag,traceSTIM)
% MAKESEGMENTS, LIKE FINDEVENTS BUT ALSO MAKES A MATRIX WITH THE SEGMENTS
% IT CAN BE USED TO FIND AND STORE UPSTATES OR POLYEVENTS: BUT MUST BE USED
% WITH DIFFERENT PARAMETERS (in both cases the accepted segments are stored
% in the upstate matrix) All data is in points and not time
%
% Upstates: for Upstate detection resonable criteria seem to be:
%minAboveTime = 1000;
%MinOnOffInt  = 2000;
%MaxIntUpstate = 250;
%SegLength     = 5000;
%preEventTime  = 1000;

% METHOD
% Finds the indexes of all events above thr.
% Calculates the diff of this vector (for continous = 1; otherise > 1)
% and uses the diff vector to get the offset-onset intervals > MaxInt.
% MaxInt is key, used to determine if 2 close thrcrossing are 1 or 2 events.
%
%%% INPUT  PARAMETERS %%%
% trace: the data
% Thr = threshold above median(trace) at which an event is defined
% MaxIntEvent = Max Int between ThrCross and still be counted as a single eventax interval within a segment
% MaxIntUpstate = Max Int between ThrCross and still be counted as a single eventax interval within a segment
% note that if MaxInt=50 and events = 1 49 98 140 it will only select event = 1;
% SegLen: length of Segment in points
% MinAboveTime: min time above thr to be accepted as a UpState. [0] to ignore.
% MinDur: min dur betweed ON and OFF be accepted as a Upstate. [0] to ignore.
% preEventTime: Amount of time to include before the onset of the detected event.
%               Also uses as posteventtime
% KurtosisFlag: 0= should be the default which ignores kurtosis selection
%               1= tries to only select UpStates based on negative kurtosis
%
%%% OUTPUT PARAMATERS s%%%
% eventtimes: vector of event onsets
% offtimes: vector of offtimes (reverse threshold crossings);
% Tabove: total time above for trace (NOT SEGMENTS)
% areaabove: total area above for trace (NOT SEGMENTS)
% SegDur (SegOff) = time in which segment trace reverse crosses Thr for the last
% Sample: matrix with all accepted segments
% DVB Modified 5/2/13: changed rest = median(trace) to = mode(round(trace))

%% SAMPLE TRACES
if (0)
    trace=zeros(1,10000); trace(10:20)=10; trace(500:510)=15; trace(520:525)=20; trace(700:800)=25; trace(1000:1005)=30; trace(1900:1910)=35; trace(2000:2300)=40; trace(3000:3005)=45; trace(3400:3405)=50;
    [eventtimes, offtimes, Tabove, areaabove, numThrX, upstate, updur, upstart, upend] = MAKESEGMENTS(trace,5,100,500,1000,300,300,100,0)
end

eventtimes=[];
potupstates=[];
offtimes=[];
upstate=[];
numThrX = 0;
updur=[];
upstart=[];
upend = [];

postTime=preEventTime;
Thr=mode(round(trace))+thr;

above =find(trace>Thr);
Tabove=length(above);
areaabove=sum(trace(above)-Thr);

if (length(above)>0)
    eventtimes = find(diff([-9999 above])>1);
    numThrX = length(eventtimes);
    eventtimes = find(diff([-9999 above])>MaxIntEvent);
    offtimes=above([eventtimes(2:length(eventtimes))-1 length(above)]);
    eventtimes=above(eventtimes);
    %Potential Upstates when Events have eventimes times > MaxIntUpstate
    potupstates = (find(diff([-9999 above])>MaxIntUpstate));
    potupstatesoff=above([potupstates(2:length(potupstates))-1 length(above)]);
    potupstates = above(potupstates);
    
end

count=0;
for j=1:length(potupstates)
    reject=0;
    start=potupstates(j);
    %Make sure Seg does not start to early or last too long
    if ( start>preEventTime && start+SegLen<length(trace) )
        seg=start-preEventTime:(start+SegLen);
        Seg=trace(seg);
        TraceStim = traceSTIM(seg);
        above =find(Seg(preEventTime:length(Seg))>Thr);
        segdur=(potupstatesoff(j)-start);
        %if (length(above)<MinAboveTime || segdur<MinOnOffInt)
        % Decide if upstate based on MinAbovetime, MinOnOffInit, and
        % eliminate if the upstate lasts longer then SegLen
        if (length(above)<MinAboveTime || segdur<MinOnOffInt || segdur>SegLen)
            reject=1;
        end
        %% REJECT BASED ON KURTOSIS
        if (KurtosisFlag)
            if (kurtosis(Seg)-3)>0  %ACCEPT NEGATIVE
                reject=1;
            end
        end
        %INCORPORATE SEGMENT INTO SAMPLE MATRIX "UPSTATE"
        if (reject==0)
            count=count+1;
            updur(count)=segdur;
            upstate(count,:)=Seg;
            upstart(count)=potupstates(j);
            upend(count) = potupstatesoff(j);
            
            %% GRAPHICS
            if (0)
                clf;
                subplot(2,1,1)
                plot(trace);
                subplot(2,1,2)
                hold on;
                plot(Seg,'r','linewidth',2)
                plot(TraceStim,'k');
                line([1 SegLen],[Thr Thr],'color','k','linestyle',':')
                line([preEventTime preEventTime],[min(trace) max(trace)],'color','k','linestyle',':','linewidth',1.5)
                line([preEventTime+segdur preEventTime+segdur],[min(trace) max(trace)],'color','k','linestyle',':','linewidth',1.5)
                str=sprintf('i=%4d, j=%4d, start/stop=%4d/%4d updur=%4d',i,j,seg(1),seg(length(seg)),segdur);
                title(str)
                drawnow;
                waitforbuttonpress
            end
            
        end
    else Seg=0; end
end
