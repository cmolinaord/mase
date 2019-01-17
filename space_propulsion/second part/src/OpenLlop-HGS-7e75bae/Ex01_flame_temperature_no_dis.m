%***********************************************************************************************************
%* HGS 1.3
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% Example 01: CH4 combustion without dissociation
% 
% Inlet: 6.25 CH4 + 25 O2 + 25.79/21 N2 at Tr
% Outlet: 12.5 H2O + 12.5 O2 + 25.79/21 N2 at Tp
 
clear; clc;
 
Tr=400;         % K
p=1;            % not rellevant 
 
% Specific enthalpy (h) of each species is obtained
hCH4=@(T) hgssingle('CH4','h',T,p);
hO2=@(T)  hgssingle('O2','h',T,p);
hN2=@(T)  hgssingle('N2','h',T,p);
hH2O=@(T) hgssingle('H2O','h',T,p);
hCO2=@(T) hgssingle('CO2','h',T,p);
 
% Combustion equation to solve: sum(hreactives) - sum(hproducts) = 0
eq=@(Tp) 6.25*hCH4(Tr)+25*hO2(Tr)+25*(79/21)*hN2(Tr)-...
         (12.5*hH2O(Tp)+12.5*hO2(Tp)+25*(79/21)*hN2(Tp)+6.25*hCO2(Tp));
 
options=optimset(...
        'Display','iter',...
        'MaxIter',4000,...
        'TolFun', 1.0e-10,...
        'TolX',1.0e-4);
 
[Tp,fval,exitflag]=fsolve(eq,1000,options);

fprintf('Adiabatic flame temperature: %.2f K \n',Tp);