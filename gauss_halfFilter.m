function [yfilt] = gauss_halfFilter(x,y,sigma)
% gaussian half filter
% for plotting firing rates from neuropixel
% https://stackoverflow.com/questions/6992213/gaussian-filter-on-a-vector-in-matlab 

    gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
    gaussFilter = gaussFilter / sum (gaussFilter); % normalize
    yfilt = conv (y, gaussFilter, 'same');
    
end

