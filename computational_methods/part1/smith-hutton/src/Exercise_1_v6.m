clc
close all
clear all

problem='DF';  %SH for Smith-Hutton, DF for diagonal flow

method='UDS';  %UDS, CDS

rho=1; 
V_0=100; %Not included inside the function so it can be easily varied for tests
beta=1; % 0 for explicit, 1 for implicit, 0.5 for Crank-Nicholson
ratiorhoGamma=[10,1e3,1e6]; %Rho/Gamma  
n=120; %Number of nodes (in x direction)
t_limit=100;
timestep=0.001;

iter_stop=1e6;
delta_conv=1e-7;


num_sol=[
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

for g=1:length(ratiorhoGamma)

    rhoGamma=ratiorhoGamma(g);
%% Mesh

[Xpos, Ypos, dx, dy, n_x, n_y]=Mesh2DGen(n, problem);        

%% Coefficients Definition

%Initial field and problem constants
alpha=10; Phi_prev=zeros(n_y,n_x); Phi_new=zeros(n_y,n_x);

%Discretization coefficients
AP=zeros(n_y,n_x); AN=zeros(n_y,n_x); AS=zeros(n_y,n_x);
AW=zeros(n_y,n_x); AE=zeros(n_y,n_x); bP=zeros(n_y,n_x);


switch problem           
    case 'SH'
[AP, AN, AS, AW, AE, bP]=BoundaryNodesSH(AP,AN,AS,AW,AE,bP,Xpos,n_x,n_y,alpha); %Boundary nodes coefficient calculations
[AP, AN, AS, AW, AE, bP]=InnerNodesSH(AP,AN,AS,AW,AE,bP,dx,dy,timestep,Xpos,Ypos,n_x,n_y,alpha,rhoGamma,rho,method); %Inner nodes coefficient calculations

    case 'DF'
[AP, AN, AS, AW, AE, bP]=BoundaryNodesDF(AP,AN,AS,AW,AE,bP,n_x,n_y); %Boundary nodes coefficient calculations
[AP, AN, AS, AW, AE, bP]=InnerNodesDF(AP,AN,AS,AW,AE,bP,dx,dy,Xpos,Ypos,n_x,n_y,rhoGamma,rho,V_0);  %Inner nodes coefficient calculations

end

%% Numerical Solution

diff=1000; iter_diff=0; n_iter=0; t=0;

while n_iter<iter_stop && diff>delta_conv && t<t_limit
    diff=0; %Difference between current and previous is reset
    iter_diff=0;
    n_iter=n_iter+1;
    
    Phi_prev=Phi_new;
    Phi_ini=Phi_prev;
    Phi_new=Phi_ini;
    
    %Recalculation of AP and bP
    switch problem
    case 'SH'
    [AP bP]=IN_Calc_SH(AP,AN,AS,AW,AE,bP,dx,dy,timestep,n_x,n_y,rho,Phi_prev,beta);
    [bP]=BN_Calc_SH(AN,AS,AW,AE,bP,n_x,n_y,Phi_prev,beta);
    case 'DF'
   [AP bP]=IN_Calc_DF(AP,AN,AS,AW,AE,bP,dx,dy,timestep,n_x,n_y,rho,Phi_prev,beta);
    [bP]=BN_Calc_DF(AN,AS,AW,AE,bP,n_x,n_y,Phi_prev,beta);
    end
    
    %Gauss-Seidel solver
    [Phi_new]=GS_solver(AP,AN,AS,AW,AE,bP,n_x,n_y,rho,Phi_ini,iter_stop,delta_conv,Phi_new,beta);
    
    for i=1:n_x
        for j=1:n_y
            iter_diff=abs(Phi_prev(j,i)-Phi_new(j,i));
            if iter_diff>diff
                diff=iter_diff;
            end     
        end
    end
    t=t+timestep;
end

%%Post-Process

figure()
    pcolor(Phi_new);
    title(strcat('Control volumes=', int2str((n_x-2)*(n_y-2)),', \rho/\Gamma=', num2str(rhoGamma,'%.2e')))
    xlabel('n_x')
    ylabel('n_y')
    shading flat;
    colormap jet;
    colorbar;
if problem=='DF'
    figure()
    contour(Phi_new(1:n_x,1:n_y));
    title(strcat('Control volumes=', int2str((n_x-2)*(n_y-2)),', \rho/\Gamma=', num2str(rhoGamma,'%.2e')))
    xlabel('n_x')
    ylabel('n_y')
    shading flat;
    colormap jet;
    colorbar;
end
    if problem=='SH'
        for h=1:(n_x/2+1)
            Phi_outlet(h)=Phi_new(1,(h-1)+(n_x/2));
        end
        X_outlet=0:1/(n_x/2):1;

    figure()
        plot(num_sol(:,1),num_sol(:,g+1),'ok',X_outlet,Phi_outlet,'b')
        xlabel('x')
        ylabel('\phi')
        legend('Exact Num. Solution',method)
        grid
        grid minor
    end
end