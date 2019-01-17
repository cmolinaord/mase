%***********************************************************************************************************
%* HGS 1.3
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% Example 05b: V2 Turbopump example (liquid H2O2 adiabatic decomposition)
%
% 80% H2O2 mass fraction
% vaporization heat of H2O2 @ 80%: 420 cal/g = 420*4.18 kcal/kg
% 
% Inlet:  molH2O2.H2O2 + molH2O.H2O
% Outlet: mollH2O2. (H2O + (1/2)O2 ) + molH2O a Tp
 
clear; clc;
 
R=8.314*1e-3;   % kJ/molK
Tr=300;         % K
pp=8;           % bar (products' pressure)

molH2O2=0.8*(1/34e-3)   % molH2O2 for kg of reactives
molH2O=0.2*(1/18e-3)    % molH2O for kg of reactives

% Computation of enthalpy, pressure is not rellevant
hH2O2=@(T) hgssingle('H2O2','h',T,1);
hO2=@(T)  hgssingle('O2','h',T,1);
hH2O=@(T) hgssingle('H2O','h',T,1);

% Equation to be solved
eq=@(Tp) molH2O2*hH2O2(Tr)+molH2O*hH2O(Tr)-420*4.18...
         -(molH2O2*(hH2O(Tp)+0.5*hO2(Tp))+molH2O*hH2O(Tp));
 
options=optimset(...
        'Display','iter',...
        'MaxIter',4000,...
        'TolFun', 1.0e-10,...
        'TolX',1.0e-4);
  
[Tp,fval,exitflag]=fsolve(eq,1000,options);
fprintf('Tp = %.2f K \n',Tp);

% Computation of the turbopump power
% first, gamma coeficient is needed
[Cp,~,MM,~,gamma,~]=hgsprop({'H2O' 'O2'},[molH2O2+molH2O,molH2O2*0.5],Tp,pp);

p2=1; % pressure at turbine outlet
nus=0.85

T2s=Tp*(p2/pp)^( (gamma-1)/gamma )

deltaT=nus*(Tp-T2s)

T2=Tp-deltaT

% Cp is in kJ/molK, and needs to be converted to kJ/kgK:
% MM in g/mol
Cp=Cp/(MM*1e-3)

w=Cp*deltaT % kJ/kg of propulsant
