function [Fe, FG] = force_vector(data, x, Tnod, L, Te, rho_eff, A, R, Re, Tn)

	%% Loads
	l_1_fs = 1.5e+4; % Lift 1 at front Spar (N/m)
	l_1_rs = 4.5e+3; % Lift 1 at rear Spar (N/m)
	l_2_fs = 5e+3;   % Lift 2 at front Spar (N/m)
	l_2_rs = 4.5e+3; % Lift 2 at rear Spar (N/m)
	d_1_fs = 1000;   % Drag 1 at front Spar (N/m)

	M_e = 1000;      % Mass of the engine (Kg)
	g0  = 9.81;      % Gravitational acceleration (m/s^2)
	W_e = M_e * g0;  % Weight of the engine (N)

	%% Lift force (Positive Z direction)
	L_f = zeros(1,data.n_el);
	% Front Spar lift
	y_1 = x(1,2);  % Coordinate y1 at Front spar
	y_2 = x(29,2); % Coordinate y2 at Front spar
	y_3 = x(53,2); % Coordinate y3 at Front spar

	for e = 28:53 % Elements at the front Spars
		y = (x(Tnod(e+1, 1), 2)) + (x(Tnod(e, 1), 2)) / 2;
		if (28 <= e && e < 42)
			L_f(e) = 0.5 * (l_1_fs + l_2_fs + (l_1_fs - l_2_fs) * cos(pi * ((y - y_1) / (y_2 - y_1)))); % Lift dist. at the front spar
		elseif (42 <= e && e < 54)
			L_f(e) = l_2_fs * cos((pi/2) * ((y - y_2) / (y_3 - y_2)));
		end
	end

	% Rear Spar lift
	y_1 = x(5,2);  % Coordinate y1 at Rear spar
	y_2 = x(32,2); % Coordinate y2 at Rear spar
	y_3 = x(54,2); % Coordinate y3 at Rear spar

	for e = 54:79 % Elements at the rear Spars
		if e == 79
			y = (x(Tnod(e, 1), 2)) + (x(Tnod(e, 1), 2)) / 2;
		else
			y = (x(Tnod(e+1, 1), 2)) + (x(Tnod(e, 1), 2)) / 2;
		end
		if (54 <= e && e < 70)
			L_f(e) = 0.5 * (l_1_rs + l_2_rs + (l_1_rs - l_2_rs) * cos(pi * ((y - y_1) / (y_2 - y_1)))); %Lift dist. at the rear spar
		elseif(70 <= e && e <= 79)
			L_f(e) = l_2_rs * cos(pi/2 * ((y - y_2) / (y_3 - y_2)));
		end
	end

	%% Weight force (Negative Z direction)
	W_f = zeros(1,data.n_el);
	% Structure weight
	for e = 1:data.n_el
		W_f(e) = rho_eff(e) * A(e) * g0;
	end

	% Engine weight (elements 7 and 8)
	W_f(7) = W_e / (2 * L(7));
	W_f(8) = W_e / (2 * L(8));

	%% Drag force (Positive X direction)
	D_f = zeros(1,data.n_el);
	% Only in Front Spar
	y_1 = x(1,2);  % Coordinate y1 at Front spar
	y_3 = x(53,2); % Coordinate y3 at Front spar

	for e = 28:53
		y = (x(Tnod(e+1, 1), 2)) + (x(Tnod(e, 1), 2)) / 2;
		D_f(e) = d_1_fs * ((1 - ((y - y_1) / (y_3 - y_1))^2));
	end

	%% Thrust force (Negative X direction)
	T_f = zeros(1,data.n_el);
	% Only in Engines (elements 7 and 8)
	T_f(7) = Te / (2 * L(7));
	T_f(8) = Te / (2 * L(8));

	%% Assembly all the forces per element
	for e = 1:data.n_el
		q_x_prima = -D_f(e) + T_f(e);
		q_y_prima = 0;
		q_z_prima = -L_f(e) + W_f(e);
		q_prima(:,e) = [q_x_prima; q_y_prima; q_z_prima];

		q_prima_barra = R(:,:,e) * q_prima(:,e);

		Fe_prima(:,e) = [(q_prima_barra(1)*L(e))/2 (q_prima_barra(2)*L(e))/2 (q_prima_barra(3)*L(e))/2 0 (-q_prima_barra(3)*(L(e)^2))/12 (q_prima_barra(2)*(L(e)^2))/12 (q_prima_barra(1)*L(e))/2 (q_prima_barra(2)*L(e))/2 (q_prima_barra(3)*L(e))/2 0 (q_prima_barra(3)*(L(e)^2))/12 (-q_prima_barra(2)*(L(e)^2))/12]';

		Fe(:,e) = Re(:,:,e)' * Fe_prima(:,e);

	end

	%% Assembly into the global force vector
	Fe = Fe';
	FG = zeros(data.ndof,1);
	for e = 1:data.n_el
		for i = 1:data.nnod * data.ngl
			I = Tn(e,i);
			FG(I) = Fe(e,i);
		end
	end

end
