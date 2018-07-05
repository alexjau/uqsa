function plot_S_i(S_i,input,parNames, lim, cMap)
% plot_S_i(S_i,parNames, lim)
% Plots an area plot based on a matrix of sensitivity indices S_i as well as a sample
% of output trajectories. The parameters are sorted based on the max value
% of the S_i for all output points.
%
% INPUT:
%      S_i: A m x p matrix containing sesnsitivity indices, where m is the number
%           of parameters and p the number of input/output points
%    input: array with p number of input points
% parNames: The parameter names in the same order as the parameters in S_i
%      lim: The limit of max(S_i) for which the parameter name is shown in
%           the legend.
%     cMap: Colormap to use
%
% For details see (Eriksson & Jauhiainen et al, 2018)
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

[nPars, ~]=size(S_i);

%Get max value of S_i for each parameter, considering all output points
maxSi=max(S_i,[],2); %Array with max(S_i) for each parameter
[~, sorted_idx]=sort(maxSi,'ascend');

%plot S_i
set(gca,'ColorOrder',cMap);
hold on
h=area(input, S_i(sorted_idx,:)');

hold on

% Get parameter names to be used in legend
leg=false(nPars,1); %Array with indicators
leg_idx=[];
for i=1:nPars
    if maxSi(sorted_idx(i))>lim
        leg(i)=1;
        leg_idx=[leg_idx sorted_idx(i)];
    else
        leg(i)=0;
    end
end

for i=1:nPars
     if leg(i)==0
         set(get(get(h(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
     end
 end
parNames   
legend(fliplr(h),flipud(parNames(leg_idx)), 'Interpreter', 'none');
%legend(fliplr(h),fliplr({parNames{leg_idx}}), 'Interpreter', 'none');

end

