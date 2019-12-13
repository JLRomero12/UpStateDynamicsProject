clearvars -EXCEPT MasterMETA
close all

cd('..\ConnectivityExperiments\Pyr-SST-Recordings')

Files = dir('Cell*.mat');
numFiles = length(Files);
[~,I] = sort([Files.datenum]);
Files = Files(I); 

% MasterMETA(1:numFiles) = struct('field1','field2','field3','field4','field5','field6');
 for file = 1:numFiles
     load(Files(file).name)
     
     MasterMETA(file) = CONN(1);
     clear META
 end