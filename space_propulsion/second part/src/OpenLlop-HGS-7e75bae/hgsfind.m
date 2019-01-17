function hgsfind( str )
%***********************************************************************************************************
%* HGS 1.3
%* By Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%Prints all the species containing str
%if str is empty, ie hgsfind(''), lists all the species
load('BurcatDB.mat');
lstr=lower(str);
for i=1:length(IndexCell)
    nom=lower(IndexCell{i,2});
    if ~isempty(strfind(nom,lstr)) | isempty(str)
        fprintf('%s\n',IndexCell{i,2});
    end
end

end

