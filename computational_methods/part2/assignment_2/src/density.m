function [rho_hat, rho_eff] = density(data, Tmat, m, V, mat)

M_w = 1000;	% Mass of the total wing (Kg)
M_s = 0;	% Mass of the Spars (Kg)
M_r = 0;	% Mass of the Ribs (Kg)
V_s = 0;	% Volume of the Spars (m^3)
V_r = 0;	% Volume of the Ribs (m^3)

for e = 1:data.n_el
	if Tmat(e) == 1
		M_s = M_s + m(e); % Total mass of the Spars
		V_s = V_s + V(e); % Total volume of the Spars
		rho(e) = mat(1,3); % Density of the Spars (kg/m^3)
	elseif Tmat(e) == 2
		M_r = M_r + m(e); % Total mass of the Ribs
		V_r = V_r + V(e); % Total volume of the Ribs
		rho(e) = mat(2,3); % Density of the Ribs (kg/m^3)
	end
end

rho_hat = (M_w - M_s - M_r) / (V_s + V_r); % Pseudo-density

rho_eff = rho + rho_hat;

end
