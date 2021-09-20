function Trial = Phase_design(start_time,lag,mid_point,binSize)
% returns time gap to test in ranksum test of FR during experiments
    if start_time == 0
        Trial=[];
    else
        end_time=start_time+lag;
        Trial=(mid_point+start_time/binSize):1:(mid_point+end_time/binSize);
        Trial=round(Trial);
    end
end

