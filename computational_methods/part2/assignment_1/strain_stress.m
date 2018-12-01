function [strain,stress] = strain_stress(x,Tnod,m,Tmat,Tdof,u)
	%Strain stress determines the strains and stresses for each bar

	for e = 1:length(Tmat)
		% Calculate the rotation matrix
		[R,l] = element_R_matrix(x, Tnod, e);

		for i = 1:6 % element displacement in global coordinates
			I = Tdof(e,i);
			u_e(i,1) = u(I,1);
		end

		u_p = R * u_e; % element displacement in local coordinates

		strain(e,1) = 1/l * [-1 1] * u_p;
		stress(e,1) = m(Tmat(e,1),1) * strain(e,1);

	end

end
