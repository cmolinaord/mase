function [Tc, frac, sp] = hydra(Tin, isliquid, x)
% This function computes hydrazine products temperature
%
% Inputs:
% Tin: inlet temperature (K)
% isliquid: 1 means liquid
% x: decomposition
% 
% Outputs:
% Tc: outlet (chamber) temperature (K)
% frac: outlet products vector (in molar fraction)
% sp: list of products

% Num of mol of each species as a function of x
nN2H4	= 1;
nNH3	= (4/3)*(1-x);
nN2		= (1/3)*(1+2*x);
nH2		= 2*x;

frac = [nNH3,nH2,nN2];
frac = frac/sum(frac);

Pc = 20; % Pressure, irrelevant for hgs enthalpy

% auxiliary functions for the enthalpies
hN2H4	= @(T) hgssingle('N2H4','h',T,Pc);
hNH3	= @(T) hgssingle('NH3','h',T,Pc);
hN2		= @(T) hgssingle('N2','h',T,Pc);
hH2		= @(T) hgssingle('H2','h',T,Pc);

if isliquid == 1
	Deltah_vap = 44.5; % kJ/mol N2H4
else
	Deltah_vap = 0;
end

DeltaH = ...
	@(Tp) nN2H4*(hN2H4(Tin)-Deltah_vap)...
	-( nNH3*hNH3(Tp) + nN2*hN2(Tp) + nH2*hH2(Tp) );

options = optimset(...
	'Display','off',...
	'MaxIter',4000,...
	'TolFun', 1.0e-10,...
	'TolX',1.0e-4);

% Solve hydrazine products temperature
Tc = fsolve(DeltaH,1000,options);

% List of species
sp = {'NH3','H2','N2'};
end
