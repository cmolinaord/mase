function Data=hgsDB(elids)
%***********************************************************************************************************
%* HGS 1.3 
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% For any issues with the code see the documentation manual.
%
% For internal use of the code.
%
%   This code is part of the HGS TOOLBOX
%   OpenLLOP, UPC-ETSEIAT 2014-2015

global BurcatDB

% Start up: Checks if BurcatDB exists on global workspace, and if it
% doesn't it loads it.
if (~exist('BurcatDB','var') || isempty(BurcatDB)), load BurcatDB; end

% Rebuild data from BurcatDB
len = length(elids);
Data = cell(len,1);

for i=1:len
    Data{i}=struct(...
        'M',BurcatDB{elids(i)}.M,...
        'CP_L',BurcatDB{elids(i)}.CP_L,...
        'CP_H',BurcatDB{elids(i)}.CP_H);
end

end