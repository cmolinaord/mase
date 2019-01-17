close all; clear all; clc;

g0      = 9.81; % Gravity acceleration [m/s^2]
Tc      = 311.11; % Chamber temperature [K]
Ae_At   = [50 80 100]; %linspace(10,100,10);

% Gas properties ordered in the following way:
% [Helium, Nitrogen, Freon-14, Freon-12, Ammonia, Hydrogen, Nitrous oxide]
% Molar masses [kg/mol]
M = [  4e-3;  28e-3;  88e-3;  121e-3;  17e-3;   2e-3;   44e-3];
% Gamma constants (specific heat ratio) [dimensionless]
g = [ 1.659;    1.4;   1.22;    1.14;   1.31;    1.4;    1.27];
% Gas constants [kJ/kg*K]
R = [2.0771; 0.2968; 0.0947; 0.06887; 0.4885; 4.1267; 0.18777];

for i=1:length(g)
    
    cs(i) = ceestar(R(i), g(i), Tc);
    
    for j=1:length(Ae_At)
        Pe_Pc(i,j)   = pressure_ratio(g(i), Ae_At(j));
        V_e(i,j)     = exit_velocity(g(i), M(i), Tc, Pe_Pc(i,j));

        F_div_m(i,j) = (V_e(i,j) + (Ae_At(j) * Pe_Pc(i,j) * cs(i)));
        Isp(i,j)     = F_div_m(i,j) / g0;
        
    end
    
end


figure()
plot(Ae_At, Pe_Pc);
xlabel('Area Ratio')
ylabel('Pressure Ratio')
title('Presure ratio vs Area ratio')

figure()
plot(Ae_At, V_e);
xlabel('Area Ratio')
ylabel('Exit velocity [m/s]')
title('V_e vs Area ratio')

figure()
plot(Ae_At, Isp);
xlabel('Area Ratio')
ylabel('Specific Impulse')
title('Isp vs Area ratio')