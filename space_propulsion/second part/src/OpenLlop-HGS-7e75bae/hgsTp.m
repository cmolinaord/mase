function [ Tp,np ] = hgsTp(species,nr,Tr,P,solver,Tstar,options)
%***********************************************************************************************************
%* HGS 1.3
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
% 
% Thermodynamical: Adiabatic combustion temperature (species dissociation, equilibrium).                         
% For any issues with the code see the documentation manual.
%
% Usage:
%       [ Tp,np ] = HGSTP(species,nr,Tr,P,solver,Tstar,options)
%
% Inputs:
%   species         -> Cell array with the species of the inlet mixture
%   nr               -> Vector for the number of mols of the inlet species
%   Tr [K]           -> Inlet temperature
%   P [bar]         -> Pressure of the chamber
%   solver         -> Select solver from fsolve/fzero to hgsfsolve
%   Tstar           -> Temperature for start solver iteration
%   options (optional) -> Options of the solver to be used. For Matlab solvers use the
%     OPTIMSET structure to generate the options structure. For custom
%     hgsfzero generate the following structure:
%
%       options = struct('x2',5000,'fchange',2,'epsx',1e-1,'epsy',1e-4,'maxite',200,'info',[]);
%
%     Where info can remain empty if nothing is to be set in screen.
%
% Output:
%   Tp              -> Products temperature (K)
%   np              -> Vector with product composition (mol)
%
% See also HGSEQ, HGSPROP, HGSSINGLE, HGSISENTROPIC, HGSFZERO
%
%   This code is part of the HGS TOOLBOX
%   OpenLLOP, UPC-ETSEIAT 2014-2015

% If info not inputed make it empty.
if ~exist('solver','var'), solver='fzero'; end
if ~exist('options','var'), options=[]; end
if ~exist('Tstar','var') || isempty(Tstar), Tstar=2500; end

[~,~,~,~,~,~,Hin,~,~]=hgsprop(species,nr,Tr,P); % Inlet enthalpy

    function DeltaH=DeltaH(Tprod)
        comp=hgseq(species,nr,Tprod,P);
        [~,~,~,~,~,~,Hout,~,~]=hgsprop(species,comp,Tprod,P);
        DeltaH=Hout-Hin;
    end

% Solving the problem
Tp = hgssolve(@DeltaH,Tstar,solver,options);

np=hgseq(species,nr,Tp,P);

end

