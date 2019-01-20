function [Phi_new] = GS_solver(a,bP,w,rho,Phi_ini,iter_lim,delta,Phi_new,beta)
	%Gauss-Seidel

	dif = 1000;
	iter = 0;
	aux = 0;

	while iter < iter_lim && dif > delta
		aux = 0;
		dif = 0;
		iter = iter+1;

		for i = 2:w.nx-1 %Inner nodes
			for j = 2:w.ny-1
				Phi_new(j,i) = (1/a.P(j,i))*(beta*(a.E(j,i)*Phi_ini(j,i+1)+a.W(j,i)*Phi_ini(j,i-1)+a.N(j,i)*Phi_ini(j+1,i)+a.S(j,i)*Phi_ini(j-1,i))+bP(j,i));
			end
		end

		for i = 2:w.nx-1 %Top and bottom, no corners
			%Top:
			Phi_new(w.ny,i) = (1/a.P(w.ny,i))*(beta*(a.E(w.ny,i)*Phi_ini(w.ny,i+1)+a.W(w.ny,i)*Phi_ini(w.ny,i-1)+a.S(w.ny,i)*Phi_ini(w.ny-1,i))+bP(w.ny,i));

			%Bottom:
			Phi_new(1,i) = (1/a.P(1,i))*(beta*(a.E(1,i)*Phi_ini(1,i+1)+a.W(1,i)*Phi_ini(1,i-1)+a.N(1,i)*Phi_ini(1+1,i))+bP(1,i));
		end

		for j = 2:w.ny-1 %Left and right, no corners
			%Left:
			Phi_new(j,1) = (1/a.P(j,1))*(beta*(a.E(j,1)*Phi_ini(j,1+1)+a.N(j,1)*Phi_ini(j+1,1)+a.S(j,1)*Phi_ini(j-1,1))+bP(j,1));

			%Right:
			Phi_new(j,w.nx) = (1/a.P(j,w.nx))*(beta*(a.W(j,w.nx)*Phi_ini(j,w.nx-1)+a.N(j,w.nx)*Phi_ini(j+1,w.nx)+a.S(j,w.nx)*Phi_ini(j-1,w.nx))+bP(j,w.nx));
		end

		%Corners:
		%Up-left
		Phi_new(w.ny,1) = (1/a.P(w.ny,1))*(beta*(a.E(w.ny,1)*Phi_ini(w.ny,1+1)+a.S(w.ny,1)*Phi_ini(w.ny-1,1))+bP(w.ny,1));
		%Up-right:
		Phi_new(w.ny,w.nx) = (1/a.P(w.ny,w.nx))*(beta*(a.W(w.ny,w.nx)*Phi_ini(w.ny,w.nx-1)+a.S(w.ny,w.nx)*Phi_ini(w.ny-1,w.nx))+bP(w.ny,w.nx));
		%Bottom-left:
		Phi_new(1,1) = (1/a.P(1,1))*(beta*(a.E(1,1)*Phi_ini(1,1+1)+a.N(1,1)*Phi_ini(1+1,1))+bP(1,1));
		%Bottom-right:
		Phi_new(1,w.nx) = (1/a.P(1,w.nx))*(beta*(a.W(1,w.nx)*Phi_ini(1,w.nx-1)+a.N(1,w.nx)*Phi_ini(1+1,w.nx))+bP(1,w.nx));

		for i = 1:w.nx
			for j = 1:w.ny
				aux = abs(Phi_new(j,i)-Phi_ini(j,i));
				if aux>dif
					dif = aux;
				end

			end
		end

		if dif > delta
			Phi_ini = Phi_new;
		end
	end
end
