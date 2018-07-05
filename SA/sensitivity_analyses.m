function [S_i, all_par_vars]=sensitivity_analyses(draws,outp, nBins)
% [S_i]=sensitivity_analyses(draws,outp, nBins) 
% Calculates the first order sensitivity index S_i according to the
% method described in (Eriksson & Jauhiainen et al, 2018) which was
% inspired from (Saltelli).
%
% INPUT:
% draws: an n x m matrix corresponding to the posterior distribution, where 
%        n is the number of draws (sample size) and m the namber of parameters
%  outp: is an m x p, matrix where m is the number of draws and p the number
%        of output points for the prediction made based on the posterior
%        distribution.
% nBins: Is the nuber of bins used in each dimension to discretize the
%        posterior distribution.
%
% OUTPUT:
%   S_i: is a m x p matrix with the sensitivity index for each parameter
%        and output point.
%
% For details see (Eriksson & Jauhiainen et al, 2018)
%
% REFERENCES:
%
% Eriksson & Jauhiainen et al., (2018) Uncertainty quantification, propagation and 
% characterization by Bayesian analysis combined with global sensitivity 
% analysis applied to dynamical intracellular pathway models", 
% https://www.biorxiv.org/content/early/2018/04/04/294884.full.pdf
%
% Saltelli, A., Tarantola, S., Campolongo, F., and Ratto, M. (2004). Sensitivity analysis
% in practice: a guide to assessing scientific models. John Wiley & Sons.
%
% Copyright (C) 2018 Olivia Eriksson (olivia@kth.se)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.

[nSamp,nPars]=size(draws);
[~,nPoints]=size(outp); 

%Calculate number of samples in each bin for each parameter and save
%bin-identities for each parameter sample point.
all_bin_counts=NaN(nPars,nBins);
all_bin_identities=NaN(nPars,nSamp);
for par_idx=1:nPars
    [bincounts,~,ind] = histcounts(draws(:,par_idx),nBins);
    all_bin_counts(par_idx,:)=bincounts(1:nBins);
    all_bin_identities(par_idx,:)=ind;
end

%Calculate mean of each bin for each parameter and inputvalue(supplemenatry eq S17)
all_bin_means=NaN(nPars,nPoints,nBins);
for par_idx=1:nPars
    for point_idx=1:nPoints
        for bin_idx=1:nBins
            bin_numb=find(all_bin_identities(par_idx,:)==bin_idx);
            if ~isempty(bin_numb)
                all_bin_means(par_idx,point_idx,bin_idx)=mean(outp(bin_numb,point_idx));
            else
                all_bin_means(par_idx,point_idx,bin_idx)=0;
                %Actually NaN, but set to zero, since this value is always multiplied with the zero probabilty further down
            end
        end
    end
end

%Calculate total mean (supplementary eq S18)
total_means=mean(outp); % Overall mean for each oupput point.

%Calculate total variance (supplementary eq S21)
total_vars=var(outp); % Overall variance for each output point.

%Calculate variance over bins with respect to parameter i and calculate S_i
%(supplementary eq S20 and S19)
all_par_vars=NaN(nPars,nPoints); %Variance for each parameter and output point.
S_i=NaN(nPars,nPoints);%sensitivity index for each parameter and output point.

for point_idx=1:nPoints
    for par_idx=1:nPars
        var_i=0; %variance of parameter i (i=par_idx) for current output point
        for bin_idx=1:nBins
            var_i=var_i+all_bin_counts(par_idx,bin_idx)*(all_bin_means(par_idx,point_idx,bin_idx)-total_means(point_idx)).^2;
        end
        var_i=var_i/nSamp; %(suopplementar eq S20)
        all_par_vars(par_idx,point_idx)=var_i;
        S_i(par_idx,point_idx)=all_par_vars(par_idx,point_idx)/total_vars(point_idx); %(supplementar eq S19)
    end
end
