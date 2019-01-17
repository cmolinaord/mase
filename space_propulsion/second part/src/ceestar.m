function [ cs ] = ceestar( R, g, Tc )
% R: [kJ/kgK]
% g: gamma [dimensionless]
% Tc: Chamber temperature [K]
% returns: cs (ceestar) [m/s]

Gamma   = sqrt(g) * ( 2 / (g + 1))^((g + 1) / (2 * (g - 1)));
cs      = sqrt(R * Tc) / (Gamma);

end