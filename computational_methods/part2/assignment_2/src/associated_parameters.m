function [L, A, I_y, I_z, J, V, m] = associated_parameters(data, x, Tnod, dat, Tmat, mat)

	for e=1:data.n_el

		% Length
		delta_x = x(Tnod(e,2),1) - x(Tnod(e,1),1);
		delta_y = x(Tnod(e,2),2) - x(Tnod(e,1),2);
		delta_z = x(Tnod(e,2),3) - x(Tnod(e,1),3);

		L(e) = sqrt(delta_x^2 + delta_y^2 + delta_z^2); % Associated length (m)

		% Area
		h = dat(e,1);
		if Tmat(e,1) == 1
			a = mat(1,4);
			b = mat(1,5);
			t = mat(1,6);

			rho = mat(1,3); % Density (kg/m^3)

		elseif Tmat(e,1) == 2
			a = mat(2,4);
			b = mat(2,5);
			t = mat(2,6);

			rho = mat(2,3); % Density (kg/m^3)
		end

		A(e) = ((h - t)*a) + (2*b*t); % Associated section area (m^2)

		% Inertia
		I_y(e) = (1/12)*a*((h - t )^3) + (1/6)*b*(t^3) + b*t*(0.5*h^2);
		I_z(e) = (1/12)*(h - t)*(a^3) + (1/6)*t*(b^3);

		% Torsion Constant
		J(e) = 1/3*(2*b*t^3 + h*a^3);

		% Volume
		V(e) = L(e)*A(e); % Associated volume (m^3)

		% Mass
		m(e) = rho*V(e); % Associated mass (kg)

	end

end
