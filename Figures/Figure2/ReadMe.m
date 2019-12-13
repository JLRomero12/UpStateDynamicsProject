% This is a read me file, not meant for running script

% All of the contents in this folders are files in order to build Figure 2.
% The raw data for the firing rate and the correlation comes from these
% files:

%The Raw examples are made with these figures:
% For PV:
% PV_SST_OptoStudy_New_19\ConnectivityExperiments\Pyr-PV-Recordings\C2_0320_.mat

% For SST:
% PV_SST_OptoStudy_New_19\ConnectivityExperiments\Pyr-SST-Recordings\S3C1_0517_.mat

% Pyr-ChR2PV-Paired-MasterMETA.mat: This file has paired data between
% pyramidal and parvalbumin neurons. This data is NOT to be confused with
% the data used for the optogenetic figures. This data comes from
% recordings done in 2018 when we had PV neurons transfected with ChR2 as
% opposed to Cheta or Halo (no double dipping up in this joint). This file
% was made from raw data in
% PV_SST_OptoStudy_New_19\ConnectivityExperiments\Pyr-PV-Recordings using
% the script UpStateAnalGen2.m in the Analysis folder. 

% SSTTdTomatoMETA.mat: This file has paired data between pyramidal and
% somatostatin neurons. The Raw data comes from this folder:
% PV_SST_OptoStudy_New_19\ConnectivityExperiments\Pyr-SST-Recordings and it
% was segmented using UpStateAnalGen2.m in the Analysis folder. This data
% is not from the data used in the Opto figures, SST neurons were
% transfected with tdTomato in order for the simultaneous recording to take
% place (as the boss says, no double dipping. This isn't the department's
% end of quarter party... ).

%%%Salvaged data (time are hard)%%%
%Since I got a bunch of recordings from slices that were not transfected
%enough I am taking the dual recordings from them and using them for figure
%2. We are only using the non-stimulated upstates from those recordings. 
%The files from those recordings are as follows: 

%1Trans-SSTCheta-PairedMETA.mat raw data comes from
%Juan_2018\PV_SST_OptoStudy_New_19\OptoExperiments\SOM-Cheta\SingleTrans
%and it was analysed with UpStateAnalGen2 and merged with merge.m in the
%Analysis folders. 
%
%1Trans-PVCheta-PairedMETA.mat raw data comes from
%Juan_2018\PV_SST_OptoStudy_New_19\OptoExperiments\PV-Cheta\SingleTrans
%and it was analysed with UpStateAnalGen2 and merged with merge.m in the
%Analysis folders. 
%
%1Trans-PVHalo-PairedMETA.mat raw data comes from
%Juan_2018\PV_SST_OptoStudy_New_19\OptoExperiments\PV-Halo\SingleTrans
%and it was analysed with UpStateAnalGen2 and merged with merge.m in the
%Analysis folders. 
%
%1Trans-SSTHalo-PairedMETA.mat raw data comes from
%Juan_2018\PV_SST_OptoStudy_New_19\OptoExperiments\SOM-Halo\SingleTrans
%and it was analysed with UpStateAnalGen2 and merged with merge.m in the
%Analysis folders. 