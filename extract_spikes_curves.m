%% Extract neuropixels data - Deborah SL 
close all
clear all
addpath(genpath('spikes/'))
addpath('Data_input/') % Data input = /camp/stp/babs/inputs/iacarusof/gaia.bianchini/deborah.schneider-luftman/Depth_decoding_neurons
load('data_Coliseum_trial.mat');
Trials=readtable('Experiment_types.csv');

%% Prepare meta data file
what=fieldnames(data);
Keep=what(1:7);
Data_interim=struct();
for ix=1:17
    for k=1:numel(Keep)
       t=Keep(k);
       Data_interim(ix).(t{1})=data(ix).(t{1});
    end
end

%%
for ix=1:17 %select which experiment (mouse x repeat)
    %%
    myDelay=data(ix).info_Coliseum_delay; %where the data are saved
    SC_IDs=myDelay.SC_IDs'; %IDs of neurons in the superior colliculus
    
    binSize=0.007;%0.005; %this determines how big the bins are - you can vary this from 0.0001 to whatever fits best   
    window=[-0.5 0.5];
    Time_line=-0.5:binSize:0.5-binSize;
    [~,Time_points]=size(Time_line);
    %
    Discard_neurons_all=zeros(1,length(SC_IDs));
   %%
    for c=1:length(SC_IDs) % for each neuron 
        % %
        st = myDelay.spikeTimes(myDelay.clu==SC_IDs(c));
        [psth, bins,~,~, spikeCounts, ba] = psthAndBA(st, myDelay.eventTimes, window, binSize);
        
        % % compute raster
        [TrialsN, inds] = sort(myDelay.trGroups); %sort the trial groups
        T=[TrialsN, inds];
        baSm = (ba./binSize); %obtain firing rate, dividing spike counts with binsize
        % % check if keeping neuron
        k=floor(Time_points/100);
        Discard_neuron=cell(10,1);
        [n1,~]=size(baSm);
        [n2,~]=size(myDelay.trGroups);
        tmp=myDelay.trGroups(1:min(n1,n2));
        for t=1:10 % for each trial type, average firing rate of 50 repeats
           % t=9 <-> visual only ||  % t=10 <-> audio only
            obj=baSm(tmp==t,:);
            %
            tmp_trace=mean(obj)';
            res_1=median(movmean(tmp_trace,k));
            res_11=median(movvar(tmp_trace,k));
            if res_11 == 0
                Discard_neuron{t}='no signal';
                if res_1 == 0
                    Discard_neuron{t}='dud';
                end
            end
        end
         % remove neurons with all duds
        Discard_neurons_all(c)=all(cellfun(@strcmp, Discard_neuron, repmat({'dud'},10,1)));
    end
    SC_IDs_updated=SC_IDs(~Discard_neurons_all);
    %%
    baSmSss=cell(length(SC_IDs_updated),1);%struct();
    tmpSs=cell(length(SC_IDs_updated),1);
    Bins=Time_line;%cell(length(SC_IDs_updated),1);
    %%
    %{
    for c=1:length(SC_IDs_updated)
        %%
        st = myDelay.spikeTimes(myDelay.clu==SC_IDs_updated(c));
        [psth, bins,~,~, spikeCounts, ba] = psthAndBA(st, myDelay.eventTimes, window, binSize);
        
        % % compute raster
        [TrialsN, inds] = sort(myDelay.trGroups); %sort the trial groups
        T=[TrialsN, inds];
        baSm = (ba./binSize); %obtain firing rate, dividing spike counts with binsize
        baSmSss{c}=baSm;
       % Bins{c}=bins;
        %
        [n1,~]=size(baSm);
        [n2,~]=size(myDelay.trGroups);
        tmp=myDelay.trGroups(1:min(n1,n2));
        tmpSs{c}=tmp;
    end
    %}
    %%
    Data_interim(ix).Data=myDelay;
    Data_interim(ix).Time.binSize=binSize;
    Data_interim(ix).Time.timeline=Bins;
    Data_interim(ix).Neurons.SC_IDs_updated=SC_IDs_updated;
end
%%
save('interim_data_FR.mat','Data_interim');