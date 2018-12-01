clear all
clc
%% Input values & constants

n = 14; % Number of nodes
nel = 41; % number of bars
nd = 3; % Number of dimensions
ndof = n * nd; % Number of degrees of freedom

g0 = 9.81; %Gravity constanst
M = 120; % Mass of the payload (Kg)
F0 = M*g0; % Force 
[x, Tnod, Tmat, Tdof, Tdn] = get_input_data();


%% External forces
F = zeros(1,ndof);
F(3) = -F0; % Force of the payload
F([18,21,24,27,30,33,36,39,42]) = F0/9; % Upper aerodynamical forces

%% Degrees of freedom 
vL = [1 2 3];
vR = 4:n;

%Material data

% Cross section of material (m2)
A = [pi*(1.5e-3)^2;                         % cable
    pi*(13.6e-2)^2-pi*(13.6e-2-6e-3)^2];    % bars
% Young modulus of materials (Pa)
E = [2e5;   % cable
    7e4];   % bars
% Density of materials (kg/m3)
rho = [1.5e3; % cable
    2.3e3];   % bars

Lc = 3.5; % Length of the cables (m)
Lbl = 3.5; % Length of the long bars (m)
Lbs = 1.25; % Length of the short bars (m)
Mpl = 120; % Payload mass (kg)
Mc = A(1)*Lc*rho(1); % Cable mass (kg)
Mbl = A(2)*Lbl*rho(2); % Long bar mass (kg)
Mbs = A(2)*Lbs*rho(2); % Long bar mass (kg)

%% Computation of element stiffness matrix   

K_el = zeros(ndof,ndof,nel);
for e=1:nel
    x(Tnod(e,1),1);
    l = element_length(x, Tnod, e);
    R = (1/l)*[x2-x1 y2-y1 z2-z1 0 0 0;
               0 0 0 x2-x1 y2-y1 z2-z1];
    K_1 = A(Tmat(e))*E(Tmat(e))/l * [1 -1; -1 1];
    K = R'*K_1*R;
    
    %store element matrix
    for r = 1:2*nd
        for s = 1:2*nd
            K_el(r,s,e) = K(r,s);
        end
    end
end

%% Global matrix

KG = zeros(ndof,ndof);
for e = 1:nel
    for i = 1:2*nd
        I = Tdof(e,i);
        for j = 1:2*nd
            J = Tdof(e,j);
            KG(I,J) = K_el(i,j,e);
        end
    end
end




%% Global system of equations

vr = 1:3; % Fixed DoF
vl = 4:42; % Free DoF

K_LL = KG(vl,vl);
K_LR = KG(vl,vr);
K_RL = KG(vr,vl);
K_RR = KG(vr,vr);

u_l = zeros(length(vl),1);
u_r = zeros(length(vr),1);

Wpl=M*g0;
Wc=Mc*g0;
Wbl=Mbl*go;
Wbs=Mbs*go;
Ws=Ms*go;
Wt=Mt*go;
Drag=Wt; %Terminal Velocity condition

F_ext=zeros(ndof,1);
F_ext(3)= -Wpl-(Wc/2)*4;
F_ext(6)= -(Wc/2)*5-(Wbs/2)*3;
F_ext(9)= -(Wc/2)*5-(Wbs/2)*2;
F_ext(12)= -(Wc/2)*5-(Wbs/2)*3;
F_ext(15)= -(Wc/2)*5-(Wbs/2)*2;
F_ext(18)= Drag/9-Ws/9-(Wc/2)-(Wbs/2)-(Wbl/2)*2;
F_ext(21)= Drag/9-Ws/9-(Wc/2)*2-(Wbs/2)*2-(Wbl/2);
F_ext(24)= Drag/9-Ws/9-(Wc/2)-(Wbs/2)-(Wbl/2)*2;
F_ext(27)= Drag/9-Ws/9-(Wc/2)*2-(Wbs/2)-(Wbl/2)*2;
F_ext(30)= Drag/9-Ws/9-(Wc/2)*4-(Wbs/2)*2-(Wbl/2)*6;
F_ext(33)= Drag/9-Ws/9-(Wc/2)*2-(Wbs/2)-(Wbl/2)*2;
F_ext(36)= Drag/9-Ws/9-(Wc/2)-(Wbs/2)-(Wbl/2)*2;
F_ext(39)= Drag/9-Ws/9-(Wc/2)*2-(Wbs/2)*2-(Wbl/2);
F_ext(42)= Drag/9-Ws/9-(Wc/2)-(Wbs/2)-(Wbl/2)*2;
F_ext_l=F_ext(vl);
F_ext_r=F_ext(vr);

u_l=K_LL\(F_ext_l-K_LR*u_r);
R=K_RR*u_r+K_RL*u_l-F_ext_r;

%Store displacements
u=zeros(ngl,1);
u(vr)=u_r;
u(vl)=u_l;

k=1;
for i=1:n
    [I,J]=find(Tdn==k);
    for j=1:nd
        disp(I,j)=u(k);
        k=k+1;
    end
end

%% Strains and stresses

for e=1:nel
    x1=x(Tnod(e,1),1);
    y1=x(Tnod(e,1),2);
    z1=x(Tnod(e,1),3);
    x2=x(Tnod(e,2),1);
    y2=x(Tnod(e,2),2);
    z2=x(Tnod(e,2),3);

    l_el=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
    
    R=(1/l_el)*[x2-x1 y2-y1 z2-z1 0 0 0; 0 0 0 x2-x1 y2-y1 z2-z1];
    
    %Obtain global element displacements
    for r=1:nnod*nd
        I=Tdof(e,r); %p, global degree of freedom
        d_g(r,1)=u(I,1);
    end
    
    %Calculate local element displacements
    d_l=R*d_g;  %Element displacements at local coordinates
    
    %Calculate strain and stresses
    epsilon(e)=(1/l_el)*[-1 1]*d_l; %Element strain
    
    if(Tmat(e)==1)
        
        sigma(e)=E1*epsilon(e); %Element stress
    
    elseif(Tmat(e)==2)
        
        sigma(e)=E2*epsilon(e);
         
    end
    
    stress_strain(e)=sigma(e)/epsilon(e);
    
end

plotDisp(x,Tnod,u',1)

figure(1)
plotStress(x,Tnod,stress_strain')


