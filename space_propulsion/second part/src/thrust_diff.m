function y = thrust_diff(x)

% This function computes the difference between the 
% theoretical Thrust value with the one computed using
% the equations of hydrazine and gas properties.
%
% For a given x (decomposition factor), it gives the error between
% theoretical thrust (1N) and the computed one. So this function
% is intended to give zero as result.


%% Input data

F = 1;				% Nominal Thrust (N)
Ae_At = 80;			% Nozzle expansion ratio by area
m_flow = 4.4e-4;	% Nominal mass flow (kg/s)

T_in = 52 + 273.15; % Inlet temperature (K) (Flash point temperature)
M = 0.0320452;	% Molar mass (kg/mol)
R0 = 8.314;		% Universal gas constant (J/mol/K)
d_e = 10e-3;	% Diameter of the exit area (m) (from engine diagram)
isliquid = 1;	% Hydrazine is liquid

Ae = pi*(d_e/2)^2;	% Exit area
At = Ae/Ae_At;		% Thoat area

Pc = 15; % Inlet pressure (bar)

%% Function implementation

[Tc, frac, sp] = hydra(T_in, isliquid, x); % Chamber temperature

[Cp,Cv,MM,Rg,gamma,a,H,G,S] = hgsprop(sp,frac,Tc,Pc); % Gamma, molar mass and R constant

M = MM*1e-3; % Molar mass (kg/mol)

[Pe_Pc] = pressure_ratio(gamma,Ae_At);

Ve = exit_velocity(gamma,M,Tc,Pe_Pc);

Pe = Pe_Pc*Pc*1e+5; % Convert to Pa

y = F - m_flow*Ve - Pe*Ae;

end