%Murrough Murnaghan
%Jack Greaves
%Syed Abdul Hadi Sharique

function dndx = odesyst(x_span,x)

global omega_c; global Ei

%Given values
cn = 300; %ms-1 Thermal velocity of neutrals
Gm = 6e22; %[m-2.s-1] Injected Partical Flux density
Gd = 1.1*Gm; %[m-2.s-1] Current density
me = 9.109384e-31; %kg electron mass
sigmaO = 3.6e-20 ; %m^2
mi = 2.2e-25; %kg
e = 1.60217662e-19; %Coulumbs

Gib = x(2)*x(1);
Geb = Gib-Gd;

Ei_prime = 2.5*Ei;

Gnb = Gm-Gib;
nn = Gnb/cn;

cbar = sqrt(8*x(3)/(pi*me));
Ri = sigmaO*cbar*(1+((2*x(3))/Ei))*exp(-Ei/x(3));
nu_ion = nn*Ri;

vex = Geb/x(2);

alpha = 1/16;
sigma_en = 27e-20;
nu_en = nn* sigma_en *sqrt(8*x(3)/(pi*me));
nu_e = (omega_c^2)/(nu_en + alpha*omega_c);

left = (5/3)*x(3)-mi*(x(1)^2);

% x(1)=0.9999*x(1);

dndx=[((5/3)*me*vex*nu_e*x(1)+nu_ion*((5/3)*x(3)+mi*x(1)*(x(1)-cn)-(x(1)/vex)*((2/3)*Ei_prime+(5/3)*x(3))))/left;
    (-(5/3)*me*x(2)*vex*nu_e-x(2)*nu_ion*(mi*(2*x(1)-cn)-((2/3)*Ei_prime+(5/3)*x(3))*(1/vex)))/left;
    (-(2/3)*mi*(x(1)^2)*me*vex*nu_e-mi*nu_ion*((2/3)*x(3)*(2*x(1)-cn)-(((mi*(x(1)^2)-x(3))/(mi*vex))*((2/3)*Ei_prime+(5/3)*x(3)))))/left;
    (-(5/3)*me*vex*nu_e*mi*(x(1)^2)-mi*nu_ion*((5/3)*x(3)*(2*x(1)-cn)-((x(1)^2)/vex)*((2/3)*Ei_prime+(5/3)*x(3))))/(left*e)];

end