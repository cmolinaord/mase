function [a bP] = InnerNodesDF(a,bP,w,Xpos,Ypos,rhoGamma,rho,v0)

	% This program is conceived as to use UDS schemes
	V_x = v0*cos(45*pi/180);
	V_y = v0*sin(45*pi/180);

	gamma = (rho/rhoGamma);

	% Calculation of the inner coefficients (in UDS)
	for i = 2:(w.nx-1)
		for j = 2:(w.ny-1)
			% Distance between consecutive nodes
			deltaXe = abs(Xpos(i+1)-Xpos(i));
			deltaXw = abs(Xpos(i)-Xpos(i-1));
			deltaYn = abs(Ypos(j+1)-Ypos(j));
			deltaYs = abs(Ypos(j)-Ypos(j-1));

			% aN coefficient
			Dn = (gamma*w.dx)/deltaYn;
			Fn = (rho*V_y)*w.dx;
			% Pen = Fn/Dn;
			Pen = 1; % To impose upwind
			if (-Fn) >= 0
				aux = -Fn;
			else
				aux = 0;
			end
			a.N(j,i) = Dn*Pen+aux;

			% aS coefficient
			Ds = (gamma*w.dx)/deltaYs;
			Fs = (rho*V_y)*w.dx;
			% Pes = Fs/Ds;
			Pes = 1; % To impose upwind
			if Fs >= 0
				aux = Fn;
			else
				aux = 0;
			end
			a.S(j,i) = Ds*Pes+aux;

			% aW coefficient
			Dw = (gamma*w.dy)/deltaXw;
			Fw = (rho*V_x)*w.dy;
			% Pew = Fw/Dw;
			Pew = 1; % To impose upwind
			if (Fw) >= 0
				aux = Fw;
			else
				aux = 0;
			end
			a.W(j,i) = Dw*Pew+aux;

			% aE coefficient
			De = (gamma*w.dy)/deltaXe;
			Fe = (rho*V_x)*w.dy;
			% Pee = Fe/De;
			Pee = 1; % To impose upwind
			if (-Fe) >= 0
				aux = -Fe;
			else
				aux = 0;
			end
			a.E(j,i) = De*Pee+aux;
		end
	end
end
