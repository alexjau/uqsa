function [ MolCaPerMolCaM ] = SSsolutionFigA(pars,Ca)
% [ MolCaPerMolCaM ] = SSsolutionFigA(pars,Ca)
% Copyright (C) 2018 
% Alexandra Jauhiainen (alexandra.jauhiainen@gmail.com)
% Olivia Eriksson
% Sara Maad Sasane
% Carolina Sartorius
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
%
% Input parameters 
% (1)  KD*CaM_Ca3*Ca
% (2)  KD*CaM_Ca2*Ca
% (3)  KD*CaM_Ca1*Ca
% (4)  KD*CaM*Ca
% (5)  KD*CaM_Ca4*PP2B
% (6)  KD*PP2B_CaM_Ca3*Ca
% (7)  KD*PP2B_CaM_Ca2*Ca
% (8)  KD*PP2B_CaM_Ca1*Ca
% (9)  KD*PP2B_CaM*Ca


% Free parameters
KD1 = pars(4); %[KD*CaM*Ca] = KD1
KD2 = pars(3); %[KD*CaM_Ca1*Ca] = KD2
KD3 = pars(2); %[KD*CaM_Ca2*Ca] = KD3
KD4 = pars(1); %[KD*CaM_Ca3*Ca] = KD4


MolCaPerMolCaM = zeros(length(Ca),1);
for ii=1:length(Ca)
    MolCaPerMolCaM(ii) = (Ca(ii)*(4*Ca(ii)^3 + 3*KD4*Ca(ii)^2 + 2*KD3*KD4*Ca(ii)+ KD2*KD3*KD4))/(Ca(ii)^4 + KD4*Ca(ii)^3 + KD3*KD4*Ca(ii)^2 + KD2*KD3*KD4*Ca(ii) + KD1*KD2*KD3*KD4);
end

end

