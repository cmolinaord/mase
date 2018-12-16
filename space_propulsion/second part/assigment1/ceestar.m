function c = ceestar(R, Tc, gamma)
Gamma = sqrt(gamma) .* (2./(gamma+1)).^((gamma+1)./(gamma-1));
c = sqrt(R*Tc)./Gamma;