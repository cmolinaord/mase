close all;
clear all;
clc;
%HAll EFFECT THRUSTER 
%ASSIGNMENT 3
%Murrough Murnaghan
%Jack Greaves
%Syed Abdul Hadi Sharique

% GIVEN
cn = 300; %[ms-1] Thermal velocity of neutrals
Gm = 6e22; %[m-2.s-1] Injected Partical Flux density
Gd =  1.1*Gm; %[m-2.s-1] Current density
B= 0.02 ; % [T] Magnectic field intensity
k= 1.380649e-23; % [j kg-1] Boltzmann constant 
me=9.109384e-31; %[kg] electron mass

%XENON PROPELLANT PROPERTIES 
sigmaO = 3.6e-20 ; %m^2
mi = 2.2e-25; %kg
V_i = 12.1; %V
e=1.60217662e-19; %Coulumbs

%THRUSTER PROPERTIES 
H = 15e-3; %[m] channel width
rm = 50e-3; %[m] mean radius
A=2*pi*rm*H; %[m^2] cross sectional area
L = 20e-3; %[m] channel length 

%INITIAL CONDITIONS @B
Gib=-0.01*Gm;
Geb=Gib-Gd;
Gnb=Gm-Gib;
nn=Gnb/cn;
ktEi=[0.1 0.15 0.2];
global Ei
Ei=e*V_i;
phi=0;

global omega_c
omega_c=e*B/me;

% ODE45 Solver
x_span=[0 0.05];
options = odeset('RelTol',10e-12,'AbsTol',10e-12);

for value = 1:3
    
    kTe=ktEi(value)*Ei;
    vi=-0.9999*sqrt((5/3)*kTe/mi);
    ne=Gib/vi;
    
    x0=[vi,ne,kTe,phi]; %INITIAL CONDITIONS
    [t{value},X{value}]=ode45(@odesyst,x_span,x0,options);
    
end
%% Plots
figure(1)
hold on
for value = 1:3
    Gi{value}=X{value}(:,1).*X{value}(:,2);
    plot (t{value},Gi{value})
    leg{value}=(['kTe/Ei =' num2str(ktEi(value))]);
end
title('Distribution of \Gamma_i')
legend(leg{1},leg{2},leg{3},'Location','northwest')
hold off

figure(2)
hold on
for value = 1:3
    plot (t{value},X{value}(:,2))
    leg{value}=(['kTe/Ei =' num2str(ktEi(value))]);
end
title ('Number of electrons')
legend(leg{1},leg{2},leg{3},'Location','northwest')
hold off

figure(3)
hold on
for value = 1:3
    plot (t{value},X{value}(:,3))
    leg{value}=(['kTe/Ei =' num2str(ktEi(value))]);
end
title ('Temperature distribution')
legend(leg{1},leg{2},leg{3},'Location','northwest')
hold off

figure(4)
hold on
for value = 1:3
    plot (t{value},X{value}(:,4))
    leg{value}=(['kTe/Ei =' num2str(ktEi(value))]);
end
title('\phi distribution')
legend(leg{1},leg{2},leg{3},'Location','southwest')
hold off

figure(5)
hold on
for value = 1:3
    M{value}=X{value}(:,1)/cn;
    plot (t{value},M{value}(:,1))
    leg{value}=(['kTe/Ei =' num2str(ktEi(value))]);
end
title('Mach number distribution')
legend(leg{1},leg{2},leg{3},'Location','southwest')
hold off

%% PART 2
InputPower=e*Gd*A;
c=sqrt((2*e*V_i)/mi);

for value = 1:3
    % Electron pressure
    X1=max(X{value}(:,1));
    X2=max(X{value}(:,2));
    X3=max(X{value}(:,3));
    mdot=X2*X1*A*mi;
    pe(value)=X2*X3;
    % Thrust
    F(value)=A*(mi*X2*X1)+A*pe(value);
    % Efficiency
    np(value)=(F(value)^2/mdot)/(2*InputPower);
end

