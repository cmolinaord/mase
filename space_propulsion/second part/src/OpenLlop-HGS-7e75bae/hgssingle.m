function out=hgssingle(specie,property,T,P)
%***********************************************************************************************************
%* HGS 1.3
%* By Arnau Miro, Pau Manent, Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
%
% Properties: calculates properties for an input element.
% For any issues with the code see the documentation manual.
%
% Usage:
%       out=HGSSINGLE(specie,property,T,P)
%
% Inputs:
%   specie     -> string with the species
%   property   -> string indicating the queried property, either 'h', 'g' or 's'
%   T [K]      -> Temperature of the mixture
%   p [bar]    -> Pressure of the mixture
%
% Output:
%   out   -> value of the queried property
%       > h [kJ/mol]
%       > g [kJ/mol]
%       > s [kJ/mol K]
%
% See also HGSEQ, HGSISENTROPIC, HGSPROP, HGSTP
%
%   This code is part of the HGS TOOLBOX
%   OpenLLOP, UPC-ETSEIAT 2014-2015

R = 8.3144621;  %[J mol^-1 K^-1]
P0 = 1;         %[bar]

property=lower(property);

id=hgsid(specie);
Data=hgsDB(id);

%Find the coeficients of each substance
coef_H=Data{1}.CP_H;
coef_L=Data{1}.CP_L;

%*******************************************************************************************************
%*(ENTHALPY) [kJ mol^-1]                                                                               *
%*******************************************************************************************************
%Take the coeficients
if((T>=1000)&&(T<=6000))
    h_coefs=[coef_H(6),coef_H(1:5)];
elseif((T>=200)&&(T<1000))
    h_coefs=[coef_L(6),coef_L(1:5)];
else
    error('hgssingle: ERROR. Temperature %.2f out of range.',T);
end

%Calculus of enthalpy
h = h_coefs(1) + sum(h_coefs(2:6).*(T.^(1:5))./(1:5));
h =h*R*1e-3; % kJ/mol

if property=='h'
    out = h;
    return
end

%*******************************************************************************************************
%*(ENTROPY) [kJ mol^-1 K^-1]                                                                           *
%*******************************************************************************************************
%Take the coeficients
if((T>=1000)&&(T<=6000))
    s_coefs=[coef_H(1:5),coef_H(7)];
elseif((T>=200)&&(T<1000))
    s_coefs=[coef_L(1:5),coef_L(7)];
else
    error('hgssingle: ERROR. Temperature %.2f out of range.',T);
end

%Evaluation of entropy
s0=R*( s_coefs(6) + s_coefs(1)*log(T) + sum(s_coefs(2:5).*(T.^(1:4))./(1:4)) );
s0=s0*1e-3;% kJ/mol;
s=s0-R*1e-3*log(P/P0);
if property=='s'
    out = s;
    return
end 
            
%*******************************************************************************************************
%(GIBBS FREE ENERGY) [kJ mol^-1]                                                                       *
%*******************************************************************************************************
g=h-T*s;

if property=='g'
    out = g;
    return
end

%*******************************************************************************************************

if (property ~='g')
    error('hgssingle: ERROR. %s unknown property',property);
end

end

