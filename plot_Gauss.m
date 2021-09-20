function sorted = plot_Gauss(obj,bins)
% make plots of mean FR after smoothing with gaussian filters 
    obj_gauss=zeros(size(obj));
    [repeats,~]=size(obj);
    %
    [~,Time_points]=size(bins);
    mid_point=floor(Time_points/2);
    Baseline_plot=1:1:mid_point-1;
    %
    smooth_1=0.015;
    smooth_2=0.02;
    %
    for jjj=1:repeats
        y=obj(jjj,:);
        y_new=gauss_halfFilter(bins,y,smooth_1);
        obj_gauss(jjj,:)=y_new;
    end
    tmp_trace_gauss=mean(obj_gauss)';
    tmp_trace_gauss=gauss_halfFilter(bins,tmp_trace_gauss,smooth_2);
    % % baseline normalisation 
    m_base=mean(tmp_trace_gauss(Baseline_plot));
    sd_base=std(tmp_trace_gauss(Baseline_plot))+0.5;
    sorted=(tmp_trace_gauss-m_base)./sd_base;
    sorted(sorted<0)=0;
    stderr=std((obj_gauss-m_base)./sd_base)./sqrt(length(obj_gauss));
    stderr=gauss_halfFilter(bins,stderr,smooth_2);
    %stderr(stderr<0)=0;
    figure(); niceBars(bins,sorted',stderr,'k',0.5); ylim([0 floor(repeats/2)]);
end

