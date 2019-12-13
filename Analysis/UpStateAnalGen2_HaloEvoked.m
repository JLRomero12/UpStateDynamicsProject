clear
close all

cd('C:\Users\BuonoLab\Desktop\Juan Luis\MainStudy\SOM-Halo')
Files = dir('S*.mat');
numFiles = length(Files);
META(1:numFiles,1) = struct();

CheckStim = figure('Name','CheckStim');

for file = 1:numFiles
    
    load(Files(file).name)
    CHUNK = [];
    
    if ~isempty(find(strcmp(CondOrder,'HALO_EVOKED'),1))
        UpStateIndex = find(strcmp(CondOrder,'HALO_EVOKED'));
        
        for cond = UpStateIndex
            
            ch1 = CH1{cond};
            ch3 = CH3{cond};
            dt = round(SampleRate{cond},1);
            
            preStim = 500/dt;
            postStim = 1000/dt;
            
            for i = 1:size(ch1,1)
                %% Finding Stimulation
                Stim        = find(ch3(i,:)>10);
                StimTimes   = find(diff([-999 Stim])>1);
                StimOnset   = Stim(StimTimes);
                StimOffset  = [Stim(StimTimes(2:end)-1) Stim(end)];
                
                chunk = zeros(length(StimOnset),preStim+postStim);
                
                %% Getting Chunk Aligned to Stimulation and saving
                for j = 1:length(StimOnset)
                    chunk(j,:) = ch1(i,(StimOnset(j)-preStim:StimOnset(j)+postStim-1));
                end
                
                CHUNK = [CHUNK; chunk];
            end
        end
        
        X = (1:size(CHUNK,2))*dt;
        
        figure(CheckStim)
        plot(repmat(X',1,size(CHUNK,1)),CHUNK','LineWidth',1.5)
        hold on
        rectangle('Position',[preStim*dt min(CHUNK(:))-10 preStim*dt 150],'FaceColor',[1 0.5 0 0.5])
        ylim([min(CHUNK(:))-5 max(CHUNK(:))+5])
        waitforbuttonpress
        clf
        
        META(file).Chunk = CHUNK;
        
    end
end

% save temp
%Notes: S1C1_1203_.mat looks like a really crappy cell.