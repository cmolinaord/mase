function [a, bP] = IN_Calc_DF(a,bP,w,timestep,rho,Phi_prev,beta)
	for i = 2:w.nx-1
		for j = 2:w.ny-1
			a.P(j,i) = (rho*(w.dx*w.dy)/timestep)+beta*(a.N(j,i)+a.S(j,i)+a.E(j,i)+a.W(j,i));
			bP(j,i) = (rho*(w.dx*w.dy)*Phi_prev(j,i)/timestep)+(1-beta)*((a.W(j,i)*(Phi_prev(j,i-1)-Phi_prev(j,i)))-(a.E(j,i)*(Phi_prev(j,i)-Phi_prev(j,i+1)))+(a.S(j,i)*(Phi_prev(j-1,i)-Phi_prev(j,i)))-(a.N(j,i)*(Phi_prev(j,i)-Phi_prev(j+1,i))));
		end
	end
end
