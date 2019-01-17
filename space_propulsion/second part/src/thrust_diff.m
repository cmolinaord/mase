function y = thrust_diff(x)
% This function computes the difference between the 
% theoretical Thrust value with the one computed using
% the equations of hydrazine and gas properties.
%
% For a given x (decomposition factor), it gives the error between
% theoretical thrust (1N) and the computed one. So this function
% is intended to give zero as result.


%% Manufacturer data
F = 1;				% Nominal Thrust (N)
Ae_At = 80;			% Nozzle expansion ratio by area
m_flow = 4.4e-4;	% Nominal mass flow (kg/s)
d_e = 10e-3;		% Diameter of the exit area (m) (from engine diagram)
T_in = 52 + 273.15; % Inlet temperature (K) (Flash point temperature)
Pc = 15;			% Inlet pressure (bar)
isliquid = 1;		% Hydrazine is liquid

%% Calculated data
Ae = pi*(d_e/2)^2;	% Exit area

%% Function implementation

% Chamber temperature and reaction products
[Tc, frac, sp] = hydra(T_in, isliquid, x); 

% Gamma, molar mass and R constant
[~,~,MM,~,gamma,~,~,~,~] = hgsprop(sp,frac,Tc,Pc);

% Molar mass to kg/mol
MM = MM*1e-3; 

% Pressure ratio
[Pe_Pc] = pressure_ratio(gamma,Ae_At);

% Exit velocity
Ve = exit_velocity(gamma,MM,Tc,Pe_Pc);

% Exit pressure in Pa
Pe = Pe_Pc*Pc*1e+5;

% Comparison of theoretical and calculated Thrust (must be zero)
y = F - m_flow*Ve - Pe*Ae;

end