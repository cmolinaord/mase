%***********************************************************************************************************
%* HGS 1.3 
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% Example 07: isentropic comparison with RPA software

function Ex07_isentropic_comparison_rpa

clear; clc;

shifting = true % conditions of shifting flow are assumed

species={'H','H2','H2O','H2O2','HO2','O','O2','OH'}

% Nozzle inlet mole fractions (RPA)
% for O/F ratio = 7.937 

ni_rpa=[ 0.0644043;...  % H
    0.1402066;...       % H2
    0.6176198;...       % H2O
    0.0000024;...       % H2O2
    0.0000367;...       % HO2
    0.0268578;...       % O
    0.0467048;...       % O2
    0.1041674];         % OH
    
Tc=3027.58  % K (RPA)
Pc=1        % bar (RPA)

% 1-Verification of properties at given Tc and Pc
disp('1-Verification of properties at given Tc and Pc');

[Cp,Cv,MM,Rg,gamma,a,H,~,S]=hgsprop(species,ni_rpa,Tc,Pc);
n=sum(ni_rpa);  % total number of mols in the mixture (1)
m=n*MM*1e-3;    % mixture total mass in kg

Pc
Tc
h=H/m   % kJ/kg
s=S/m   % kJ/kgK
cp=Cp/m % kJ/kgK
cv=Cv/m
Rg
MM
gamma
rho=Pc/(Rg*Tc)
a

% Results are more or less ok except cp and cv which are quite different.

% 2-Verification of the equilibrium composition given T and P
disp('Verification of the equilibrium composition given T and P');

ni_calc=hgseq(species,ni_rpa,Tc,Pc)

error_tantpercent=100*(ni_calc-ni_rpa) / norm(ni_rpa)

% 3-Verification of the insentropic expansion at a certain pressure Pt
disp('3-Verification of the insentropic expansion at a certain pressure Pt');

Pt=0.1 

% Solver options
options = struct('x2',5000,'fchange',2,'epsx',1e-1,'epsy',1e-4,'maxite',200,'info',1);

[ Tt,nt ] = hgsisentropic(species,ni_rpa,Tc,Pc,Pt,'shifting','hgsfzero',300,options)

%{
    function DeltaS=DeltaS(T)
        if shifting
            nt=hgseq(species,ni_rpa,T,Pt); % Shifting, equilibrium conditions at T,Pe
        else
        nt=ni_rpa; % frozen
        end
        [~,~,MM2,~,~,~,~,~,S2]=hgsprop(species,nt,T,Pt); % enthropy at T,Pt
        n2=sum(nt); % mixture total number of mols (1)
        m2=n2*MM2*1e-3; % mixture mass in kg
        s2=S2/m2; % kJ/kgK
        DeltaS=s2-s;
    end

Tt=fzero(@DeltaS,3000,optimset('Display','iter'))
if shifting
    nt=hgseq(species,ni_rpa,Tt,Pt) 
else
    nt=ni_rpa
end
%}

[~,~,MM2,~,~,a2,H2,~,S2]=hgsprop(species,nt,Tt,Pt);
n2=sum(nt); % mitxture total number of mols (1)
m2=n2*MM2*1e-3; % mixture mass in kg
s2=S2/m2; % kJ/kgK
fprintf('Specific entropy at nozzle inlet=%f outlet=%f kJ/kgK \n',s,s2);

h2=H2/m2;

vt=sqrt(2*1000*(h-h2)); % Enthalpy must be en J/kg !
fprintf('Nozzle outlet velocity=%f (m/s) \n',vt);
fprintf('Nozzle outlet Mach number=%f \n',vt/a2);
Is=vt/9.81; % Is (optimal expansion, Pe=Pambient)
fprintf('Specific impulse = %f (s) \n',Is);
end

