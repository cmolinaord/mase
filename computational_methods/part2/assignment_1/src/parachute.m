%  Parachute structural analysis
%  Carlos Molina 2018
%  Computational Methods - Master's Degree in Space and Aeronautical Engineering (MASE)
%  Universitat Politecnica de Catalunya (UPC)

close all
clear
clc
%% Input values & constants

n = 14; % Number of nodes
nel = 41; % number of bars
nd = 3; % Number of dimensions
ndof = n * nd; % Number of degrees of freedom

g0 = 9.81; % Gravity constanst (m/s2)
rho_air = 1.225; % Air density (kg/m3)

m_pl = 120; % Mass of the payload (Kg)
[x, Tnod, Tmat, Tdof, Tdn] = get_input_data();

% Material data
% =======================
% Upper surface material
 rhos = 1500; % kg/m3
 S = 17.5; % m2
 Ts = 0.001; % m
 Cd = 1.25;
 m_s = rhos*S*Ts; % Kg

% Material properties
radius = [0.75e-3; 6.8e-3]

% Young modulus (Pa) - Yield strength (Pa) -  Section area (m^2) - Density (kg/m3)
m = [
200000e6  300e6     pi*(radius(1)^2)*pi			1500;   % Cables
70000e6   240e6     pi*((radius(2)^2)-(5.3e-3^2))	2300];  % Bars

m_bars = zeros(nel, 1);
for e = 1:nel
	[R,len] = element_R_matrix(x, Tnod, e);
	% Element mass = density * section area * length
	m_bars(e,1) = m(Tmat(e,1),4) * m(Tmat(e,1),3) * len;
end
% Total mass
m_tot = sum(m_bars) + m_s + m_pl;

% Global stiffnes computation
K = global_stiffness(x,Tnod,m,Tmat,Tdof);

% Global_force: Drag + Weight
% Weight
W = m_tot*g0;

% Terminal velocity reached when Drag = Weight
Vt = sqrt(2*W/(rho_air*S*Cd));

% Drag
D = 0.5*rho_air*(Vt^2)*Cd*S;

% Forces per node applied in each DoF
f = zeros(ndof,1);
for i = 1:n
	% Weight per node due the connected bars
	m1 = 0.5 * sum(m_bars(Tnod(:,1) == i));
	m2 = 0.5 * sum(m_bars(Tnod(:,2) == i));
	f(i*3,1) = (m1 + m2)*g0; % (N)
	if i > 5
		% Surface weight addded and Drag substracted
		f(i*3,1) = f(i*3,1) + (m_s/9)*g0 - (D/9);
	end
end
f(3,1) = f(3,1) + m_pl*g0; % Payload weight added

% Fixed DoF
vr = [1 2 3 28 29 20];

[u,r] = global_displacements_reactions(K,f,vr);

%% Visualization
% ======================
[strain,stress] = strain_stress(x,Tnod,m,Tmat,Tdof,u);
plotParachute([1],u,strain,stress,x,Tnod,Tmat)

%% Other results
% =======================

% Stress safety factor
SF = m(Tmat,2)./abs(stress);
SF_min = min(SF);
SF_min_idx = find(SF == SF_min);

% Critical buckling stress
m(Tmat,2)
% Critical stress (Pa)
Stress_cr = zeros(nel, 1);
for e = 1:nel
	[R,len] = element_R_matrix(x, Tnod, e);
	% Element mass = Young Mod
	Stress_cr(e) = m(Tmat(e),1) * pi^2 / (len/0.5*radius(Tmat(e)))^2;
end
% Only compute maximum buckling stress for rigid bars (Infinite for cables)
Stress_cr(Tmat==1) = Inf;
SF_buckling = Stress_cr ./ abs(stress);
SF_buckling_min = min(SF_buckling);
SF_buckling_min_idx = find(SF_buckling==SF_buckling_min);

printf('The worst buckling Safety ratio is archieved in bar %i with a value of %1.1e\n', SF_buckling_min_idx, SF_buckling_min)
