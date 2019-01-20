function [bP] = BN_Calc_SH(a,bP,w,Phi_prev,beta);

	%% Update of the bottom outlet
	for i = (w.nx/2):w.nx
		bP(1,i) = (1-beta)*(a.N(1,i)*Phi_prev(2,i));
	end

	%% Update of the corners

	%Up-left:
	bP(w.ny,1) = (1-beta)*(a.S(w.ny,1)*Phi_prev(w.ny-1,1)+a.E(w.ny,1)*Phi_prev(w.ny,2));

	%Up-right:
	bP(w.ny,w.nx) = (1-beta)*(a.S(w.ny,w.nx)*Phi_prev(w.ny-1,w.nx)+a.W(w.ny,w.nx)*Phi_prev(w.ny,w.nx-1));

	%Bottom-left:
	bP(1,1) = (1-beta)*(a.N(1,1)*Phi_prev(2,1)+a.E(1,1)*Phi_prev(1,2));

	%Bottom-right:
	bP(1,w.nx) = (1-beta)*(a.N(1,w.nx)*Phi_prev(2,w.nx)+a.W(1,w.nx)*Phi_prev(1,w.nx-1));

end
