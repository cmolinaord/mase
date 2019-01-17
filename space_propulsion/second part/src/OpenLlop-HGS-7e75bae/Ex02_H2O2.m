%***********************************************************************************************************
%* HGS 1.3
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% Example 02: Adiabatic H2O2 decomposition
%
% Inlet: H2O2 
% Outlet: H2O + (1/2)O2 at Tp
 
clear; clc;
 
R=8.314*1e-3;   % kJ/molK
Tr=300;         % K
p=1;            % not rellevant
 
% Enthalpy of each species 
hH2O2=@(T) hgssingle('H2O2','h',T,1);
hO2=@(T)  hgssingle('O2','h',T,1);
hH2O=@(T) hgssingle('H2O','h',T,1);

% Equation to be solved
eq=@(Tp) hH2O2(Tr)-...
         (hH2O(Tp)+0.5*hO2(Tp));
 
options=optimset(...
        'Display','iter',...
        'MaxIter',4000,...
        'TolFun', 1.0e-10,...
        'TolX',1.0e-4);
 
 
[Tp,fval,exitflag]=fsolve(eq,1000,options);

fprintf('Dissociation temperature: %.2f K \n',Tp);