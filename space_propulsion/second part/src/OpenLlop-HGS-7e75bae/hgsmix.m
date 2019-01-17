function Out=hgsmix(Data,property,T,P,n)
%***********************************************************************************************************
%* HGS 1.3
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% For any issues with the code see the documentation manual
%
% Calculates enthalpy, enthropy and Gibbs free energy for a input element
% vector with its elements properties stored in Data.
%
%   This code is part of the HGS TOOLBOX
%   OpenLLOP, UPC-ETSEIAT 2014-2015

%***********************************************************************************************************
%*1 CONSTANTS                                                                                              *
%***********************************************************************************************************

R = 8.3144621;  %[J mol^-1 K^-1]
P0 = 1;         %[bar]

% Preallocating
len_data = length(Data);
coef_H = zeros(len_data,7);
coef_L = zeros(len_data,7);
H = zeros(len_data,1);
S0 = zeros(len_data,1);
S = zeros(len_data,1);
G = zeros(len_data,1);
CP = zeros(len_data,1);

for i=1:len_data

%***********************************************************************************************************
%*3 INPUTS                                                                                                 *
%***********************************************************************************************************
property=lower(property);


%Find the coeficients of each substance
coef_H(i,:)=Data{i}.CP_H;
coef_L(i,:)=Data{i}.CP_L;
       
%***********************************************************************************************************
%*4 FUNCTIONS                                                                                              *
%***********************************************************************************************************

    %*******************************************************************************************************
    %*4.1 h (ENTHALPY) [kJ mol^-1]                                                                         *
    %*******************************************************************************************************
        
       %Take the coeficients
        if((T>=1000)&&(T<=6000))
            h_coefs=[coef_H(i,6),coef_H(i,1:5)];
        elseif((T>=200)&&(T<1000))
            h_coefs=[coef_L(i,6),coef_L(i,1:5)];
        else
            error('hgsmix:ERROR. Temperature %.2f out of range.',T);
        end

        %Calculus of enthalpy
        h = h_coefs(1) + sum(h_coefs(2:6).*(T.^(1:5))./(1:5));
        H(i)=h*R*1e-3; % kJ/mol

    %*******************************************************************************************************
    %*4.2 s (ENTROPY) [kJ mol^-1]                                                                *
    %*******************************************************************************************************
        if (n(i) > 0)
            Pi=P*n(i)/sum(n);
            %Take the coeficients
            if((T>=1000)&&(T<=6000))
                s_coefs=[coef_H(i,1:5),coef_H(i,7)];
            elseif((T>=200)&&(T<1000))
                s_coefs=[coef_L(i,1:5),coef_L(i,7)];
            else
                error('hgsmix:ERROR. Temperature %.2f out of range.',T);
            end

            %Calculus of entropy
            s0=R*( s_coefs(6) + s_coefs(1)*log(T) + sum(s_coefs(2:5).*(T.^(1:4))./(1:4)) );
            S0(i)=s0*1e-3;% kJ/mol
            s=s0-R*log(Pi/P0);
            S(i)=s*1e-3; % kJ/mol
        else
            S(i)=0; % kJ/mol
            S0(i)=0;
        end
            
    %*******************************************************************************************************
    %*4.3 g (GIBBS FREE ENERGY) [kJ mol^-1]                                                                         *
    %*******************************************************************************************************
        G(i)=H(i)-T*S(i);
        
    %*******************************************************************************************************
    %*4.4 cp (CP) [J mol^-1 K^-1]                                                                         *
    %*******************************************************************************************************
        %Take the coeficients
        if((T>=1000)&&(T<=6000))
            cp_coefs=coef_H(i,1:5);
        elseif((T>=200)&&(T<1000))
            cp_coefs=coef_L(i,1:5);
        else
            error('hgsmix:ERROR. Temperature %.2f out of range.',T);
        end

        %Calculus of Cp
        cp = R*( cp_coefs(1) + sum(cp_coefs(2:5).*(T.^(1:4))) );
        CP(i)=cp*1e-3; % kJ/(mol*K) 
end

switch property
    case 'h'
        Out=H;
    case 's'
        Out=S;
    case 's0'
        Out=S0;
    case 'g'
        Out=G;
    case 'cp'
        Out=dot(n,CP)/sum(n);
end

end

