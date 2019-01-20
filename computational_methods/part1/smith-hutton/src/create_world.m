function [posX, posY, w] = create_world(n, problem)

	switch problem
		case 'SH'
		w.nx = n; %Horizontal axis number of nodes
		w.ny = w.nx/2; %vertical axis number of nodes
		x_ini = -1; %Initial X value for mesh
		x_final = 1; %Final X value for mesh
		y_ini = 0; %Initial Y value for mesh
		y_final = 1; %Final Y value for mesh
		case 'DF'
		w.nx = n;
		w.ny = w.nx;
		x_ini = 0;
		x_final = 1;
		y_ini = 0;
		y_final = 1;
		otherwise
		error ('Wrong problem')
	end

	%Per each axis, two frontier nodes
	w.dx = (x_final-x_ini)/(w.nx-2);
	w.dy = (y_final-y_ini)/(w.ny-2);

	%Nodes position is defined
	posX = zeros(w.nx,1);
	posX(1) = x_ini;
	posX(w.nx) = x_final;
	for i = 2:(w.nx-1)
		posX(i) = (posX(1)+w.dx/2)+w.dx*(i-2);
	end

	posY = zeros(w.ny,1);
	posY(1) = y_ini;
	posY(w.ny) = y_final;
	for i = 2:(w.ny-1)
		posY(i) = (posY(1)+w.dy/2)+w.dy*(i-2);
	end

end
