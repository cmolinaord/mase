function [uint, fint] = local_DisRotFor(data, Tn, F, u, Re, Ke)

	for e = 1:data.n_el
		for i = 1:data.nnod * data.ngl
			I = Tn(e,i);
			fe(e,i) = F(I);
			ue(e,i) = u(I);
		end
	end
	ue = ue';

	for e = 1:data.n_el
		fint(:,e) = Re(:,:,e) * Ke(e) * ue(:,e);
		uint(:,e) = Re(:,:,e) * ue(:,e);
	end

end
