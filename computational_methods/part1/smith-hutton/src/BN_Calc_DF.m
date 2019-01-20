function [bP] = BN_Calc_DF(a,bP,w,Phi_prev,beta)

	bP(w.ny,1) = (1-beta)*(a.S(w.ny,1)*Phi_prev(w.ny-1,1)+a.E(w.ny,1)*Phi_prev(w.ny,2));  %Upper-left corner
	bP(w.ny,w.nx) = (1-beta)*(a.S(w.ny,w.nx)*Phi_prev(w.ny-1,w.nx)+a.W(w.ny,w.nx)*Phi_prev(w.ny,w.nx-1));  %Upper-right corner
	bP(1,1) = (1-beta)*(a.N(1,1)*Phi_prev(2,1)+a.E(1,1)*Phi_prev(1,2)); %Lower-left corner
	bP(1,w.nx) = (1-beta)*(a.N(1,w.nx)*Phi_prev(2,w.nx)+a.W(1,w.nx)*Phi_prev(1,w.nx-1));  %Lower-right corner

	for i = 1:w.nx
		bP(w.ny,i) = (1-beta)*(a.S(1,i)*Phi_prev(w.ny-1,i)); %Upper outlet
		bP(i,w.nx) = (1-beta)*(a.W(i,w.nx)*Phi_prev(i,w.nx-1)); %Right outlet
	end

end
