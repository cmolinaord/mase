function index=hgsid(element)
%***********************************************************************************************************
%* HGS 1.3
%* By Arnau Miro, Pau Manent, Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
%* HGS. Expansion: Returns index of input element regarding hgs database.
%
% For any issues with the code see the documentation manual.
%
% For internal code usage.
%
%   This code is part of the HGS TOOLBOX
%   OpenLLOP, UPC-ETSEIAT 2014-2015

global IndexCell

% Start up: Checks if BurcatDB exists on global workspace, and if it
% doesn't it loads it.
if (~exist('IndexCell','var') || isempty(IndexCell)), load BurcatDB; end

% Search for index in the IndexCell
search = ismember(IndexCell(:,2),element);

% If not found return an error
if(all(search==0)), error('ERROR: Unknown %s',element); end % MANEL

index = IndexCell{search,1};

end
