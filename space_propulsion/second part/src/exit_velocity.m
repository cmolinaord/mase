function [ V_e ] = exit_velocity( g, M, Tc, Pe_Pc )
% ve: nozzle exit velocity
% g: gamma [dimensionless]
% M: Molar masses [kg/mol]
% Tc: Chamber temperature [K]
% P0, Pe: chamber pressure, outlet pressure
% returns: ve (nozzle exit velocity) [m/s]

R0 = 8.314; % Universal gas constant [J/molK]

V_e = sqrt(((2 * g * R0 * Tc) / ((g - 1) * M)) * ...
    (1 - Pe_Pc^((g - 1) / g)));

end