function [ T2,n2,v2,M2 ] = hgsisentropic( species,n1,T1,P1,P2,eql,solver,Tstar,options )
%***********************************************************************************************************
%* HGS 1.3
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% Expansion: Isentropic expansion (frozen or shifting flow)
% For any issues with the code see the documentation manual.
%
% Usage:
%       [T2,n2,v2,M2]=HGSISENTROPIC( species,n1,T1,P1,P2,eql,solver,Tstar,options )
%
% Inputs:
%   species       -> Cell array with the species of the inlet mixture
%   n1                -> Vector for the number of mols of each inlet species
%   T1 [K]          -> Inlet temperature
%   P1 [bar]      -> Pressure of the chamber
%   P2 [bar]      -> Outlet pressure
%   eq1              -> 'shifting' or 'frozen'
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
%   T2              -> Products temperature (K)
%   n2              -> Vector for the number of mols of each outlet species
%   v2              -> outlet velocity m/s, assuming inlet velocity is 0 
%   M2             -> outlet Mach number, assuming inlet velocity is 0
%
% See also HGSEQ, HGSPROP, HGSSINGLE, HGSTP, HGSFZERO
%
%   This code is part of the HGS TOOLBOX
%   OpenLLOP, UPC-ETSEIAT 2014-2015

% If info not inputed make it empty.
if ~exist('solver','var'), solver='hgsfzero'; end
if ~exist('options','var'), options=[]; end
if ~exist('Tstar','var') || isempty(Tstar), Tstar=T1/2; end

if strcmpi(eql,'shifting')==0 && strcmpi(eql,'frozen')==0
    error('eql should be ''shifting'' or ''frozen''');
end

[~,~,MM1,~,~,~,H1,~,S1]=hgsprop(species,n1,T1,P1); % Inlet properties

m1=sum(n1)*MM1*1e-3;
h1=H1/m1;

% Solving the problem
T2 = hgssolve(@DeltaS,Tstar,solver,options);

if strcmpi(eql,'shifting')==1
    n2=hgseq(species,n1,T2,P2);
else
    n2=n1;
end


[~,~,MM2,~,~,a2,H2,~,S2]=hgsprop(species,n2,T2,P2); % Inlet properties

m2=sum(n2)*MM2*1e-3;
h2=H2/m2;

v2=sqrt(2*1000*(h1-h2)); % Enthalpy must be en J/kg !

M2=v2/a2;

function DeltaS=DeltaS(T)
    if strcmpi(eql,'shifting')==1
        n2=n1;
        n2=hgseq(species,n1,T,P2); % Shifting 
    else
        n2=n1; % frozen
    end
    [~,~,~,~,~,~,~,~,S2]=hgsprop(species,n2,T,P2); 
    DeltaS=S2-S1;
end

end

