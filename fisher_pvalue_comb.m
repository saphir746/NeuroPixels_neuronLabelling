function [h_fisher,group_pval] = fisher_pvalue_comb(pvals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to compute Fisher's p-value for meta-analysis
%    stolen online
%       https://uk.mathworks.com/matlabcentral/fileexchange/65327-function-to-compute-fisher-s-p-value-for-meta-analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pvals is a vector which holds p-values
    chi_vals = -2.*log(pvals);
    group_pval = 1 - chi2cdf(sum(chi_vals),2*length(pvals));
    h_fisher =(group_pval<0.05);
end