% Computational Fluid Dynamics
% Convection-Diffusion Equation (Smidth-Hutton)
% MASE (Master's Degree in Space and Aeronautical Engineering)
% Carlos Molina
% January 2019

clc
clear
close all

problem = 'DF';  %SH for Smith-Hutton, DF for diagonal flow
method = 'CDS';  %UDS, CDS

rho = 1;
V_0 = 100; %Not included inside the function so it can be easily varied for tests
beta = 0.5; % 0 for explicit, 1 for implicit, 0.5 for Crank-Nicholson
ratiorhoGamma = [10,1e3,1e6]; % Rho/Gamma
n = 50; % Number of nodes (in x direction)
t_limit = 100;
timestep = 0.001;

N_prog = 40; % Number of points for the progress bar

iter_stop = 1e4;
delta_conv = 1e-5;

num_sol = [
0.0 1.989 2.0000 2.000;
0.1 1.402 1.9990 2.000;
0.2 1.146 1.9997 2.000;
0.3 0.946 1.9850 1.999;
0.4 0.775 1.8410 1.964;
0.5 0.621 0.9510 1.000;
0.6 0.480 0.1540 0.036;
0.7 0.349 0.0010 0.001;
0.8 0.227 0.0000 0.000;
0.9 0.111 0.0000 0.000;
1.0 0.000 0.0000 0.000];

tic;
% Creating struct for results
res.time = zeros(1,length(ratiorhoGamma));
res.iter = zeros(1,length(ratiorhoGamma));
res.error = zeros(1,length(ratiorhoGamma));

% Creating title results line
fprintf('Prob\tMthd\tStep\tRho/Gamma\tProgress bar')
for k = 1:3+N_prog-12; fprintf(' '); end
fprintf('\tTime\tIters\tFreq\terr\n')

for g = 1:length(ratiorhoGamma)
	fprintf('%2s\t%3s\t%2i\t%1.1e\t\t[', problem, method, g, ratiorhoGamma(g))
	rhoGamma = ratiorhoGamma(g);

	%% Mesh
	[Xpos, Ypos, w] = create_world(n, problem);

	%% Coefficients Definition

	%Initial field and problem constants
	alpha = 10;
	Phi_prev = zeros(w.ny,w.nx);
	Phi_new = zeros(w.ny,w.nx);

	%Discretization coefficients
	a.P = zeros(w.ny,w.nx); a.N = zeros(w.ny,w.nx); a.S = zeros(w.ny,w.nx);
	a.W = zeros(w.ny,w.nx); a.E = zeros(w.ny,w.nx); bP = zeros(w.ny,w.nx);

	switch problem
		case 'SH'
		[a, bP] = BoundaryNodesSH(a,bP,Xpos,w,alpha); %Boundary nodes coefficient calculations
		[a, bP] = InnerNodesSH(a,bP,w,timestep,Xpos,Ypos,alpha,rhoGamma,rho,method); %Inner nodes coefficient calculations

		case 'DF'
		[a, bP] = BoundaryNodesDF(a,bP,w); %Boundary nodes coefficient calculations
		[a, bP] = InnerNodesDF(a,bP,w,Xpos,Ypos,rhoGamma,rho,V_0);  %Inner nodes coefficient calculations

	end

	%% Numerical Solution

	error = 1000; iter_error = 0; iter = 0; t = 0;

	while iter < iter_stop && error > delta_conv && t < t_limit
		error = 0; %Difference between current and previous is reset
		iter_error = 0;
		iter = iter+1;

		Phi_prev = Phi_new;
		Phi_ini = Phi_prev;
		Phi_new = Phi_ini;

		%Recalculation of AP and bP
		switch problem
			case 'SH'
			[a, bP] = IN_Calc_SH(a,bP,w,timestep,rho,Phi_prev,beta);
			[bP] = BN_Calc_SH(a,bP,w,Phi_prev,beta);
			case 'DF'
			[a, bP] = IN_Calc_DF(a,bP,w,timestep,rho,Phi_prev,beta);
			[bP] = BN_Calc_DF(a,bP,w,Phi_prev,beta);
		end

		%Gauss-Seidel solver
		[Phi_new] = GS_solver(a,bP,w,rho,Phi_ini,iter_stop,delta_conv,Phi_new,beta);

		for i = 1:w.nx
			for j = 1:w.ny
				iter_error = abs(Phi_prev(j,i)-Phi_new(j,i));
				if iter_error > error
					error = iter_error;
				end
			end
		end
		t = t+timestep;
		% Progress bar
		if mod(iter,iter_stop/N_prog) == 0; fprintf('#'); end
	end
	for k = 1:N_prog-floor(round(iter/iter_stop*N_prog,1))
	fprintf('_')
	end
	fprintf(']\t')

	%%Post-Process

	figure()
	pcolor(Phi_new);
	title(strcat('Control volumes = ', int2str((w.nx-2)*(w.ny-2)),', \rho/\Gamma = ', num2str(rhoGamma,'%.2e')))
	xlabel('nx')
	ylabel('ny')
	shading flat;
	colormap jet;
	colorbar;
	axis('equal')
	
	if problem == 'DF'
		figure()
		contour(Phi_new(1:w.nx,1:w.ny));
		title(strcat('Control volumes = ', int2str((w.nx-2)*(w.ny-2)),', \rho/\Gamma = ', num2str(rhoGamma,'%.2e')))
		xlabel('nx')
		ylabel('ny')
		shading flat;
		colormap jet;
		colorbar;
	end
	if problem == 'SH'
		for h = 1:(w.nx/2+1)
			Phi_outlet(h) = Phi_new(1,(h-1)+(w.nx/2));
		end
		X_outlet = 0:1/(w.nx/2):1;

		figure()
		plot(num_sol(:,1),num_sol(:,g+1),'ok',X_outlet,Phi_outlet,'b')
		title(strcat('Control volumes = ', int2str((w.nx-2)*(w.ny-2)),', \rho/\Gamma = ', num2str(rhoGamma,'%.2e')))
		xlabel('x')
		ylabel('\phi')
		legend('Exact Num. Solution',method)
		grid
		grid minor
	end

	% Time computation
	time_ = toc;
	if g == 1
		res.time = toc;
	else
		res.time(g) = toc-res.time(g-1);
	end

	% Results saving
	res.error(g)	= error;
	res.iter(g)		= iter;
	fprintf('%1.2f\t', res.time(g))
	fprintf('%i\t%3.0f\t%1.1e\n', iter, iter/res.time(g), error)
end
