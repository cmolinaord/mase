function Kg = global_stiffness(x,Tnod,m,Tmat,Tdof)
	%GLOBAL STIFFNESS Returns the global stiffness matrix of the system

	% Calculate the element stiffness
	for e = 1:length(Tmat)
		[R, l] = element_R_matrix(x, Tnod, e);

		K = (m(Tmat(e,1),3)*m(Tmat(e,1),1)/l) * [1 -1; -1 1];
		Kp = R'*K*R;

		for r = 1:6
			for s = 1:6
				Ke(r,s,e) = Kp(r,s);
			end
		end
	end

	% Calculate the global stiffness
	Kg = zeros(42,42);
	for e = 1:length(Tmat)
		for i = 1:6
			I = Tdof(e,i);
			for j = 1:6
				J = Tdof(e,j);
				Kg(I,J) = Kg(I,J) + Ke(i,j,e);
			end
		end
	end

end
