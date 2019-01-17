%***********************************************************************************************************
%* HGS 1.1 
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% Example 06: LOX - LH2 combustion
%
% Temperature of stoichiometric combustion of H2 O2
% reactives inlet as satured liquit at 10 bar
%
% O2 (NIST) hv(404.36 K)-hl(119.62 K)=14.3753 kJ/mol
% H2 (NIST) hv(413.96K)-hl(31.39K)=10.9495 kJ/mol

function [Tp,np]=Ex06_LOX_combustion

clear; clc;

species={'H2','O2' , 'H2O','H','O','OH'};
nr=[2;1;0;0;0;0]; % mol

% Enthalpy of liquid O2 at Tsat 10 bar (kJ/mol)
hO2=hgssingle('O2','h',404.36,10) -14.3753; 

% Enthalpy of liquid H2 at Tsat 10 bar (kJ/mol)
hH2=hgssingle('H2','h',413.96,10) -10.9495; 

% Enthalpy of stoichiometric mixture
HinLIQ=2*hH2+1*hO2;

P=10
 
    function DeltaH=DeltaH(Tprod)
        comp=hgseq(species,nr,Tprod,P);
        [~,~,~,~,~,~,Hout,~,~]=hgsprop(species,comp,Tprod,P);
        DeltaH=Hout-HinLIQ;
    end
    
Tp=fzero(@DeltaH,3000,optimset('Display','iter'));
np=hgseq(species,nr,Tp,P);

end