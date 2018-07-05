function [ MolCaMPerMolPP2B ] = SSsolutionFigB(pars, totalCaM,totalPP2B)
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
KD9 = pars(5); %[KD*CaM_Ca4*PP2B] = KD9

% Parameters defined by rules
KD8 = (KD4 * KD9)/pars(6); %[KD*CaM_Ca3*PP2B] = KD8
KD7 = (KD3 * KD8)/pars(7); %[KD*CaM_Ca2*PP2B] = KD7
KD6 = (KD2 * KD7)/pars(8); %[KD*CaM_Ca1*PP2B] = KD6
KD5 = (KD1 * KD6)/pars(9); %[KD*CaM*PP2B] = KD5

MolCaMPerMolPP2B = zeros(length(totalCaM),1);
for ii=1:length(totalCaM)
    MolCaMPerMolPP2B(ii) = -(totalCaM(ii)*(KD5 - totalPP2B + totalCaM(ii) - (KD5^2 ...
    + 2*KD5*totalPP2B + 2*KD5*totalCaM(ii) + totalPP2B^2 ...
    - 2*totalPP2B*totalCaM(ii) + totalCaM(ii)^2)^(1/2)))/(totalPP2B*(KD5 ...
    + totalPP2B - totalCaM(ii) + (KD5^2 + 2*KD5*totalPP2B + 2*KD5*totalCaM(ii) ...
    + totalPP2B^2 - 2*totalPP2B*totalCaM(ii) + totalCaM(ii)^2)^(1/2)));
end

end