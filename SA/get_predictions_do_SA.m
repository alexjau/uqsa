%SCRIPT to get predictions based on the posterior sample and do sensitivity analysis. 
%The prediction corresponds to phenotype 4, which is used as a test case
%
% For details and references concerning the sensitifity analysis see the function:
% [S_i]=sensitivity_analyses(draws,outp, nBins) 
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


%Load posterior distribution
post=load('Draws-Phenotype123-Scale1000.mat');
[nSamples,nPars]=size(post.draws);

%input=Ca
nPoints=26;
nBins=20;
Ca=linspace(2,4,nPoints);

%output/prediction
outp=NaN(nSamples, nPoints);

%As a test case use phenotype 4 as prediction: 
for i=1:nSamples
outp(i,:)=SSsolutionFigC(10.^post.draws(i,:), 10.^Ca, 300, 3);
end


%Do sensitivity analysis on the prediction:
[S_i, all_par_vars]=sensitivity_analyses(post.draws,outp, nBins);

%plot result
subplot(2,1,1)
plot(Ca,outp(1:10:end,:));
subplot(2,1,2)
plot_S_i(S_i,Ca,post.parNames,0.0,cat(1,cool(54),jet(12)))
