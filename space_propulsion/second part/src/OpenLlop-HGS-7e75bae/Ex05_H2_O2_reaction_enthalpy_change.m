%***********************************************************************************************************
%* HGS 1.3
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% Example 05: Adiabatic H2 / O2 reaction
%
% Inlet: H2, O2 
% Outlet: H2O + (1/2)O2 at Tp

clear
close 

species={'H2','O2','H2O','H','O','OH'};
Tr=350; % K 
P=10; % bar
nr=[2;1;0;0;0;0]; % mol
 
[~,~,~,~,~,~,Hin,~,~]=hgsprop(species,nr,Tr,P)
 
T=linspace(300,5000,10);
for i=1:length(T) 
    comp=hgseq(species,nr,T(i),P);
    [~,~,~,~,~,~,Hout(i),~,~]=hgsprop(species,comp,T(i),P);    
end
 
plot(T,Hout,'-or',T,Hin*ones(length(T),1),'-ob')