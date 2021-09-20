close all
clear all
addpath(genpath('spikes/'))
addpath('Data_input/')
load('interim_data_FR.mat');
Trials=readtable('Experiment_types.csv');
%%
what=fieldnames(Data_interim);
Keep=what(1:7);
Data_new=struct();
for ix=1:17
    for k=1:numel(Keep)
        t=Keep(k);
       Data_new(ix).(t{1})=Data_interim(ix).(t{1});
    end
end
%%
for ix=1:17 %select which experiment (mouse x repeat)
    % %
    TMP=Data_interim(ix);
    %
    Data=TMP.Data;
    %
    bins=TMP.Time.timeline;
    binSize=TMP.Time.binSize;
    window=[-0.5 0.5];
    SC_IDs_updated=Data.SC_IDs;%TMP.Neurons.SC_IDs_updated;
    % % define important time lines
    [~,Time_points]=size(bins);
    mid_point=floor(Time_points/2);
    Baseline=round(0.3/binSize):1:mid_point-1;
    % %
    Sorted=zeros(length(SC_IDs_updated),Time_points,10);
    %Stderr=zeros(length(SC_IDs_updated),Time_points,10);
    Labels=zeros(length(SC_IDs_updated),10);
    % %
    for c=1:length(SC_IDs_updated)
        % % compute raster
        st = Data.spikeTimes(Data.clu==SC_IDs_updated(c));
        [psth, bins,~,~, spikeCounts, ba] = psthAndBA(st, Data.eventTimes, window, binSize);
        [TrialsN, inds] = sort(Data.trGroups); %sort the trial groups
        T=[TrialsN, inds];
        baSm = (ba./binSize); %obtain firing rate, dividing spike counts with binsize
        [n1,~]=size(baSm);
        [n2,~]=size(Data.trGroups);
        tmp=Data.trGroups(1:min(n1,n2));
        % %
       % baSm=Data{c};
       % tmp=Tmps_all{1};
        % %
        labels=zeros(1,10);
        sorted_sub=zeros(Time_points,10);
       % stderr_sub=zeros(Time_points,10);
        % %
        for t=1:10 % for each trial type, average firing rate of 50 repeats
           % % t=1 <-> audio-visual simul. || t=9 <-> visual only ||  % t=10 <-> audio only
            obj=baSm(tmp==t,:);
            [repeats,~]=size(obj);
            %{
            tmp_trace=mean(obj)';
            tmp_std=std(obj)./sqrt(length(obj));
            % sanity check plot
            figure(); niceBars(bins,tmp_trace',tmp_std,'k',0.5);
            %}
            % % gaussian filtering for plottingPha  
            line_plot=plot_Gauss(obj,bins);
            close all
            sorted_sub(:,t)=line_plot;
            %stderr_sub(:,t)=stderr;
            %%
            % for all tests - this is relevant
            X=obj(:,Baseline);
            lag=0.075;
            Trial_0=mid_point:1:(mid_point+lag/binSize);
            switch t
                case {1, 9 ,10}
                    start_time=0;
                case 2
                    start_time=0.025;
                case 3
                    start_time=0.05;
                case 4
                    start_time=0.075;
                case 5
                    start_time=0.1;
                case 6
                    start_time=0.15;
                case 7
                    start_time=0.2;
                case 8
                    start_time=0.3;
            end
            Trial_1 = Phase_design(start_time,lag,mid_point,binSize);
            Trial=[Trial_0, Trial_1];
            Trial=unique(Trial);
            Y=obj(:,Trial);
            %
            PPs=zeros(1,repeats);
            for jjj=1:repeats
                PPs(jjj)= ranksum(X(jjj,:),Y(jjj,:));%signrank(X(jjj,:),Y(jjj,:));%
            end
            PPs=PPs(~isnan(PPs));%PPs(isnan(PPs))=1;%
            [~,Pp1]=fisher_pvalue_comb(PPs);
           % if Pp1 < 0.05/10
            labels(t)=Pp1;%1;%Trials{t,2}{1};
            %end
            %
        end
        % %
         Sorted(c,:,:)=sorted_sub;
        % Stderr(c,:,:)=stderr_sub;
         Labels(c,:)=fdr_bh(labels,0.05,'pdep')';
    end
    % %
    cgo = clustergram(Labels,'Cluster','column','Colormap',redbluecmap); %''Standardize','Row'
    cgo.ColumnPDist = 'correlation';
    set(cgo,'Linkage','complete','Dendrogram',3)
    set(cgo,'ColumnLabels',Trials{:,2},'rowLabels',SC_IDs_updated);
    Perc_detected=sum(sum(Labels,2)>0)/length(SC_IDs_updated)*100;
    Perc_detected=round(Perc_detected,1);
    mytitle='Round '+string(ix)+', mouse '+string(TMP.MouseID)+'\n('+string(Perc_detected)+'% fired)';
    title = addTitle(cgo,compose(mytitle),'Color','red');
    addXLabel(cgo,'Experiments');
    addYLabel(cgo,'Neuron ID');
    %
    plot(cgo);
    set(gcf,'position',[200 200 600 800])
    saveas(gcf,'Cluster_plots/clustergram_'+string(ix)+'.png')
    % %
    Data_new(ix).line=Sorted;
   % Data_new(ix).Stdev=Stderr;
    Data_new(ix).Neurons.SC_IDs=SC_IDs_updated;
    Data_new(ix).Neurons.Labels=Labels;
end
%%
close all
% %
save('reformatted_data_MFR.mat','Data_new');
