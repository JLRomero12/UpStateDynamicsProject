function [SimilMatrix] = UpstateCorr(UpState, UpDur, dt, preeventpoints, GRAPHICS)
% Function to do the pairwise correlation between Upstates segemented from MAKESEGMENTS

postEventTime = 500; %100 ms
postEventPoints = postEventTime/dt;


[numSamples numSamplePoints]=size(UpState);
SimilMatrix = zeros(numSamples,numSamples);

if (numSamples>1)
    for i=1:numSamples-1
        for j=i+1:numSamples
            clf; hold on;
            maxLen=max(UpDur(i),UpDur(j)) + postEventTime + preeventpoints;
            if (maxLen>numSamplePoints)
                maxLen=numSamplePoints;
            end
            window=1:maxLen;
            trace1=UpState(i,window );
            trace2=UpState(j,window );
            corr = corrcoef(trace1,trace2); corr=corr(1,2);
            if (GRAPHICS)
                plot(trace1,'b'); plot(trace2,'r'); set(gca,'xlim',[0 numSamplePoints])
                dummy=sprintf('[%3d %3d] corr=%6.3f',i,j,corr);
                title(dummy)
                drawnow
                %waitforbuttonpress
            end
            SimilMatrix(i,j)=corr; SimilMatrix(j,i)=corr;
        end
    end
end

SimilMatrix(1:numSamples+1:end) = 1;

