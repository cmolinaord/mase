clear all
clc

%% Input data

% Nodal coordinates
x = [%   X       Y       Z
     0.000,  0.000,  0.000; % 1
    -0.725, -0.731,  3.250; % 2
     0.725, -0.731,  3.250; % 3
     0.725,  0.731,  3.250; % 4
    -0.725,  0.731,  3.250; % 5
     2.900, -1.450,  5.630; % 6
     2.900,  0.000,  5.820; % 7
     2.900,  1.450,  5.630; % 8
     0.000,  1.450,  6.340; % 9
     0.000,  0.000,  6.530; % 10
     0.000, -1.450,  6.340; % 11
    -2.900,  1.450,  5.630; % 12
    -2.900,  0.000,  5.820; % 13
    -2.900, -1.450,  5.630; % 14
];

% Bar connectivities
Tnod = [%   A     B
         1     2  % 1
         1     3  % 2
         1     4  % 3
         1     5  % 4
         2     3  % 5
         3     4  % 6 
         4     5  % 7 
         5     2  % 8
         2     4  % 9
         2    10  % 10
         2    11  % 11
         2    13  % 12
         2    14  % 13
         3     6  % 14
         3     7  % 15
         3    10  % 16
         3    11  % 17
         4     7  % 18
         4     8  % 19
         4     9  % 20
         4    10  % 21
         5     9  % 22
         5    10  % 23
         5    12  % 24
         5    13  % 25
        14    11  % 26
        11    10  % 27
        14    10  % 28
        14    13  % 29
        13    10  % 30
        13    12  % 31
        12    10  % 32
        12     9  % 33
         9    10  % 34
         9     8  % 35
         8     7  % 36
        10     7  % 37
        10     8  % 38
        11     6  % 39
         6     7  % 40
        10     6  % 41
];

% Material connectivities
Tmat = [% Mat. index
         1  % 1
         1  % 2
         1  % 3
         1  % 4
         2  % 5
         2  % 6 
         2  % 7 
         2  % 8
         2  % 9
         1  % 10
         1  % 11
         1  % 12
         1  % 13
         1  % 14
         1  % 15
         1  % 16
         1  % 17
         1  % 18
         1  % 19
         1  % 20
         1  % 21
         1  % 22
         1  % 23
         1  % 24
         1  % 25
         2  % 26
         2  % 27
         2  % 28
         2  % 29
         2  % 30
         2  % 31
         2  % 32
         2  % 33
         2  % 34
         2  % 35
         2  % 36
         2  % 37
         2  % 38
         2  % 39
         2  % 40
         2  % 41
];

% Bar connectivities in terms of deg of freedom
Tdof = [%    A            B
         1   2   3   4    5   6  % 1
         1   2   3   7    8   9  % 2
         1   2   3   10   11  12  % 3
         1   2   3   13   14  15  % 4
         4   5   6   7    8   9  % 5
         7   8   9   10   11  12  % 6 
         10  11  12  13   14  15  % 7 
         13  14  15  4    5   6  % 8
         4   5   6   10   11  12  % 9
         4   5   6   28   29  30  % 10
         4   5   6   31   32  33  % 11
         4   5   6   37   38  39  % 12
         4   5   6   40   41  42  % 13
         7   8   9   16   17  18  % 14
         7   8   9   19   20  21  % 15
         7   8   9   28   29  30  % 16
         7   8   9   31   32  33  % 17
         10  11  12  19   20  21  % 18
         10  11  12  22   23  24  % 19
         10  11  12  25   26  27  % 20
         10  11  12  28   29  30  % 21
         13  14  15  25   26  27  % 22
         13  14  15  28   29  30  % 23
         13  14  15  34   35  36  % 24
         13  14  15  37   38  39  % 25
         40  41  42  31   32  33  % 26
         31  32  33  28   29  30  % 27
         40  41  42  28   29  30  % 28
         40  41  42  37   38  39  % 29
         37  38  39  28   29  30  % 30
         37  38  39  34   35  36  % 31
         34  35  36  28   29  30  % 32
         34  35  36  25   26  27  % 33
         25  26  27  28   29  30  % 34
         25  26  27  22   23  24  % 35
         22  23  24  19   20  21  % 36
         28  29  30  19   20  21  % 37
         28  29  30  22   23  24  % 38
         31  32  33  16   17  18  % 39
         16  17  18  19   20  21  % 40
         28  29  30  16   17  18  % 41
];

%Material data

A1=pi*(1.5e-3)^2; %Cross section of the cables [m^2]
E1=2e+5; %Young Modulus of the cables [Pa]
rho1=1.5e+3; %Density of the cables [kg/m^3]
A2=pi*(13.6e-2)^2-pi*(13.6e-2-6e-3)^2; %Cross section of the cables [m^2]
E2=7e+4; %Young Modulus of the bars [Pa]
rho2=2.3e+3; %Density of the bars [kg/m^3]
Mpl=120; %Payload mass [kg]

%Structure data

nd=3;       %Problem dimension
nel=41;     %Total number of bars
n=14;       %Total number of nodes
ngl=n*nd;   %Total degrees of freedom
nnod=2;     %Number of nodes in a bar

%% Computation of element stiffness matrix   
    
for e=1:nel
    x1=x(Tnod(e,1),1);
    y1=x(Tnod(e,1),2);
    z1=x(Tnod(e,1),3);
    x2=x(Tnod(e,2),1);
    y2=x(Tnod(e,2),2);
    z2=x(Tnod(e,2),3);
    
    l_el=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
    
    Q=(1/l_el)*[x2-x1 y2-y1 z2-z1 0 0 0; 0 0 0 x2-x1 y2-y1 z2-z1];
    
    if(Tmat(e)==1)
        
        K_1=((A1*E1)/l_el)*[1 -1; -1 1];
    
    elseif(Tmat(e)==2)
        
         K_1=((A2*E2)/l_el)*[1 -1; -1 1];
         
    end
    
    K=Q'*K_1*Q;
    
    %store element matrix
    for r=1:nnod*nd
        for s=1:nnod*nd
            K_el(r,s,e)=K(r,s);
        end
    end
end

%% Global matrix

KG=zeros(ngl,ngl);

for e=1:nel
    for i=1:nnod*nd
        I=Tdof(e,i);
        for j=1:nnod*nd
            J=Tdof(e,j);
            KG(I,J)=KG(I,J)+K_el(i,j,e);
        end
    end
end

%% Global system of equations

vr=[1 2 3];
vl=4:42;

K_LL=KG(vl,vl);
K_LR=KG(vl,vr);
K_RL=KG(vr,vl);
K_RR=KG(vr,vr);


