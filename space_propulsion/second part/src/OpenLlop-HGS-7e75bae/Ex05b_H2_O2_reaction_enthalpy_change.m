%***********************************************************************************************************
%* HGS 1.3
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% Example 05b: Adiabatic H2 / O2 reaction using hgsTp
%
% Inlet: H2, O2 
% Outlet: H2O + (1/2)O2 at Tp

clear

species={'H2','O2','H2O','H','O','OH'};
Tr=350; % K 
P=10; % bar
nr=[2;1;0;0;0;0]; % mol

% Solver options
options = optimset('Algorithm','interior-point','Display','iter','TolFun',1e-8,'Tolx',1e-8);

[Tp,np]=hgsTp(species,nr,298,1,'fzero',[],options)
np/sum(np)