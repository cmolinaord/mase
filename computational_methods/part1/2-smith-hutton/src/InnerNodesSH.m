function [a bP] = InnerNodesSH(a,bP,w,timestep,Xpos,Ypos,alpha,rhoGamma,rho,method)

	gamma = (rho/rhoGamma);

	%% Velocity field
	[V_x_e,V_x_w,V_y_n,V_y_s] = SH_vel(Xpos,Ypos,w);

	%Calculation of the inner coefficients
	for i = 2:(w.nx-1)
		for j = 2:(w.ny-1)
			%Distance between consecutive nodes
			deltaXe = abs(Xpos(i+1)-Xpos(i));
			deltaXw = abs(Xpos(i)-Xpos(i-1));
			deltaYn = abs(Ypos(j+1)-Ypos(j));
			deltaYs = abs(Ypos(j)-Ypos(j-1));

			Dn = (gamma*w.dx)/deltaYn;
			Fn = (rho*V_y_n(j,i))*w.dx;
			Ds = (gamma*w.dx)/deltaYs;
			Fs = (rho*V_y_s(j,i))*w.dx;
			Dw = (gamma*w.dy)/deltaXw;
			Fw = (rho*V_x_w(j,i))*w.dy;
			De = (gamma*w.dy)/deltaXe;
			Fe = (rho*V_x_e(j,i))*w.dy;

			Pe_n = Fn/Dn;
			Pe_s = Fs/Ds;
			Pe_w = Fw/Dw;
			Pe_e = Fe/De;

			switch method
				case 'CDS'
				A_Pe_n = 1-0.5.*abs(Pe_n);
				A_Pe_s = 1-0.5.*abs(Pe_s);
				A_Pe_w = 1-0.5.*abs(Pe_w);
				A_Pe_e = 1-0.5.*abs(Pe_e);
				case 'UDS'
				A_Pe_n = 1;
				A_Pe_s = 1;
				A_Pe_w = 1;
				A_Pe_e = 1;
				%             case 'HDS'
				%         a(1) = 0;
				%         a(2) = 1-0.5.*abs(Pen);
				%         A_Pe_n = max(a);
				%      case 'PLDS'
				%         a(1) = 0;
				%         a(2) = (1-0.1.*abs(Pe)).^5;
				%         A_Pe = max(a);
				%      case 'EDS'
				%         if abs(Pe) =  = 0
				%             A_Pe = 1;
				%         else
				%         A_Pe = abs(Pe)./(exp(abs(Pe))-1);
				%         end
				otherwise
				error('Wrong metohd')
			end


			%aN coefficient

			%Pen = Fn/Dn;
			Pen = A_Pe_n; %To impose upwind
			if (-Fn) >= 0
				aux = -Fn;
			else
				aux = 0;
			end
			a.N(j,i) = Dn*Pen+aux;

			%aS coefficient
			Ds = (gamma*w.dx)/deltaYs;
			Fs = (rho*V_y_s(j,i))*w.dx;
			%Pes = Fs/Ds;
			Pes = A_Pe_s; %To impose upwind
			if Fs >= 0
				aux = Fn;
			else
				aux = 0;
			end
			a.S(j,i) = Ds*Pes+aux;

			%aW coefficient
			Dw = (gamma*w.dy)/deltaXw;
			Fw = (rho*V_x_w(j,i))*w.dy;
			%Pew = Fw/Dw;
			Pew = A_Pe_w; %To impose upwind
			if (Fw) >= 0
				aux = Fw;
			else
				aux = 0;
			end
			a.W(j,i) = Dw*Pew+aux;

			%aE coefficient
			De = (gamma*w.dy)/deltaXe;
			Fe = (rho*V_x_e(j,i))*w.dy;
			%Pee = Fe/De;
			Pee = A_Pe_e;
			if (-Fe) >= 0
				aux = -Fe;
			else
				aux = 0;
			end
			a.E(j,i) = De*Pee+aux;

		end
	end

	return
