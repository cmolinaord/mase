function [ ve ] = exit_velocity(gamma, M, Tc, Pe_P0)
    
% ve: nozzle exit velocity
% g: gamma
% M: kg/mol
% T0: K
% P0, Pe: chamber pressure, outlet pressure

 R0 = 8.314; %universal gas constant [J/molK]
 
 ve = sqrt(2*gamma*R0*Tc/((gamma-1)*M)*(1-(Pe_P0)^((gamma-1)/gamma)));

end