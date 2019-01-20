function [a, bP] = IN_Calc_SH(a,bP,w,timestep,rho,Phi_prev,beta);

	source = 0; %a.s per problem specifications

	%The inner nodes are updated

	for i = 2:w.nx-1
		for j = 2:w.ny-1
			%aP coefficient update
			a.P(j,i) = (rho*(w.dx*w.dy)/timestep)+beta*(a.N(j,i)+a.S(j,i)+a.E(j,i)+a.W(j,i))-source*(w.dx+w.dy);

			%bP coefficient update
			bP(j,i) = (rho*(w.dx*w.dy)*Phi_prev(j,i)/timestep)+(source*(w.dx+w.dy))+(1-beta)*((a.W(j,i)*(Phi_prev(j,i-1)-Phi_prev(j,i)))-(a.E(j,i)*(Phi_prev(j,i)-Phi_prev(j,i+1)))+(a.S(j,i)*(Phi_prev(j-1,i)-Phi_prev(j,i)))-(a.N(j,i)*(Phi_prev(j,i)-Phi_prev(j+1,i))));
		end
	end
end
