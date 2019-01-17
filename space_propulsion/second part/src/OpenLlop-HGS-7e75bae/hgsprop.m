function [Cp,Cv,MM,Rg,gamma,a,H,G,S]=hgsprop(species,n,T,P)
%***********************************************************************************************************
%* HGS 1.3
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% Properties: properties of a mixture of gases.
% For any issues with the code see the documentation manual.
%
% Usage:
%       [Cp,Cv,MM,Rg,gamma,a,H,G,S]=HGSPROP(species,n,T,P)
%
% Inputs:
%   species         -> Cell array with the species of the mixture
%   n               -> Vector for the number of mols of the species
%   T [K]           -> Temperature of the mixture
%   p [bar]         -> Pressure of the mixture
%
% Output:
%   Cp [kJ/molK]  	-> Heat coef at constant pressure
%   Cv [kJ/molK]    -> Heat coef at constant volume
%   MM [g/mol]      -> Molar mass of the mixture
%   Rg [kJ/kgK]     -> Gas constant of the mixture
%   gamma           -> Gamma coef of the mixture
%   a  [m/s]        -> Sound speed of the mixture
%   H  [kJ]         -> Enthalpy of the mixture
%   G  [kJ]         -> Gibbs free energy of the mixture
%   S  [kJ/K]         -> Entropy of the mixture
%
% See also HGSEQ, HGSISENTROPIC, HGSSINGLE, HGSTP
%
%   This code is part of the HGS TOOLBOX
%   OpenLLOP, UPC-ETSEIAT 2014-2015

R = 8.3144621/1000;  %[kJ mol^-1 K^-1]

% Search the identificators and the data of each specie
len_species = length(species);
id_species = zeros(len_species,1);
for i=1:len_species
    id_species(i)=hgsid(species{i});
end

Data=hgsDB(id_species);

Cp=hgsmix(Data,'cp',T,P,n);% kJ/(mol*K)
Cv=Cp-R;% kJ/(mol*K)
gamma=Cp/Cv;
M = zeros(len_species,1);
for i=1:len_species
    M(i)=Data{i}.M;
end
MM=dot(n,M)/sum(n); %[g/mol]
Rg=1000*R/MM; %[kJ/(KgK)]
a=sqrt(gamma*1e3*Rg*T);

H=dot(n,hgsmix(Data,'h',T,P,n));% kJ
G=dot(n,hgsmix(Data,'g',T,P,n));% kJ
S=dot(n,hgsmix(Data,'s',T,P,n));% kJ

end
