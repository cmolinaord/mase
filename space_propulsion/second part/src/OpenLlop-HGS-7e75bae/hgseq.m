function [neq,deltaG] = hgseq(species,n0,T,p,options)
%***********************************************************************************************************
%* HGS 1.3
%* By Arnau Miro, Pau Manent and Manel Soria
%
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% Equilibrium: Chemical equilibrium of a group of species.                                       
% For any issues with the code see the documentation manual.
%
% Usage:
%       [neq,deltaG] = HGSEQ(species,n0,T,p,options)
%
% Inputs:
%   species -> Cell array with the species of the reaction
%               eg. {'H2' 'O2' 'H2O' 'H' 'O' 'OH'}
%   n   	-> Vector for the number of mols of the reactives
%               eg. [2 1 0 0 0 0 ] means 2 mol of H2, one of O2 and zero of
%               the rest
%   T [K]   -> Temperature in which the reaction occurs
%   p [bar] -> Pressure in which the reaction occurs
%
%  Options: Options for the solver are given in the form of optimset
%  function. For more information, refer to FMINCON and OPTIMSET.
%
% Output:
%   neq   	-> Number of mols at equilibrium for given species
%               order of index for given example:
%               H2  O2  H2O  H  O  OH
%               1   2   3    4  5  6
%  deltaG -> Gibbs free energy of the equilibrium mixture
%
% See also HGSISENTROPIC, HGSPROP, HGSSINGLE, HGSTP, FMINCON, OPTIMSET
%
%   This code is part of the HGS TOOLBOX
%   OpenLLOP, UPC-ETSEIAT 2014-2015

% If no options are provided, use a default
if ( ~exist('options','var') || isempty(options) )
    options=optimset('Algorithm','interior-point','Display','off','TolFun',1e-9,'Tolx',1e-9);
end

%Take the identification number of every specie ang keep in the id vector
len = length(species);
id = zeros(len,1);
for i=1:len
    id(i)=hgsid(species{i});
end
%Take the data stored in the BurcatDB MAT-file for the species introduced
Data=hgsDB(id);

%Function to minimize, free energy of Gibbs
gtotal=@(n) dot(n,hgsmix(Data,'g',T,p,n));

%Disequality restrictions Ax<b.
%The concentration of each specie must to be 0 or a positive number, necer
%negative
A=eye(len,len);
A=-A;
b=zeros(len,1);

%Equality restriction Ax=b
%The total mols of each element must to conserve
%system_matrix
[Aeq,beq]=system_matrix(species,n0);

%Solution to the problem
[n,deltaG,exitflag]=fmincon(gtotal,n0,A,b,Aeq,beq,[],[],[],options);
neq=n;

if exitflag~=1 && exitflag~=2 % Manel
    error('hgseq: fmincon failed to find equilibrium')
end

end

function [M,b]=system_matrix(species,n_o)
%***********************************************************************************************************
%* HGS 1.1 
%* By Arnau Miro, Pau Manent, Eva Astrid Jara
%* Supervised by Manel Soria and Ramon Carreras
%* LLOP, ETSEIAT UPC          
%***********************************************************************************************************
%
% For any issues with the code see the documentation manual
%
% For internal use of the code
%
%NOMENCLATURE
%el -> element
%el2-> element2, in the string, position after element
%el3-> element3, in the string, position after element2
%su -> substance, specie
%f  -> row of the system matrix
%j  -> roll of the system matrix
%m  -> multiplicator, feed the matrix
%elecomp -> String vector that contain the elements apeared before
%example, elecomp={'H', 'O', 'Al'}, it no contents species, only elements

elecomp={' '};
B=zeros(10,length(species));
c=zeros(1,10);

for j=1:length(species)
    %Read every specie
    su=cell2mat(species(j));
    
    % Treat some special cases (Arnau)
    switch (su)
        case '(Jet-A(g))'
            su = 'C12H23';
        case '(Jet-A(L))'
            su = 'C12H23';
        case '(air)' % Full composition
            su = 'O2N2Ar';
    end
        
    %Search every element in the specie
    for i=1:length(su)
        el=su(i);
        % The content of the parentheses is not important in the matrix building.
        % Stop iterating and pass to next element.
        if (double(el)==40), break; end

        % If 'el' doesn't begin by capital letter, it is not a chemical element.
        % Move on to the next character
        if (double(el)<65 && double(el)>90), continue; end

        % Beginning by capital letter, four cases are available: 
        %Case 1: The element is in the last position of the script like O in H2O
        if (i==length(su))                                              
            m=1;
        else
        %Case 2: The next letter is a capital letter, the element is only a letter, like H in HO
            el2=su(i+1);
            if (double(el2)>=65 && double(el2)<=90)                     
                m=1;
        %Case 3: The next letter is a number, The element has only a capital letter 
        % but is multiplied by the number that follows, like 2 in H2O
            elseif (double(el2)>=48 && double(el2)<=57)                 
                m=double(el2)-48;
                if (i+1<length(su))
                    el3=su(i+2);
                    %Case 3.2: The next letter is a number, The element has only a 
                    % capital letter but is multiplied by the number that follows, like 2 in H2O
                    if (double(el3)>=48 && double(el3)<=57)                 
                        m2=double(el3)-48;
                        m=10*m+m2;
                    end
                end
        %Case 4: The next letter is a lowercase letter. 
        % The element has two letters like Br in Br2
            elseif(double(el2)>=97 && double(el2)<=122)
                el=[el, el2];
                m=1;
                %Is this element in the last position of the string? 
                %If the answer is yes, m=1, otherwise we need to check 
                %if the next letter is a number of a Capital letter
                if (i+1<length(su))                  
                    el3=su(i+2);
                    %The next letter is a number, change the multiplier m by this number
                    if (double(el3)>48 && double(el3)<57)               
                        m=double(el3)-48;
                        if (i+2<length(su))
                            el4=su(i+3);
                            if (double(el4)>=48 && double(el4)<=57)         
                                m2=double(el4)-48;
                                m=10*m+m2;
                            end
                        end
                    end
                end
            end
        end
        % Now we have an element and it's quantity. We build the system
        % matrix
        % Rows are different elements
        % Rolls are different species
        
        % First of all we search the corresponding row of the element.
        if (double(el(1))>=65 && double(el(1))<=90)
            aux={el};
            trobat=0;
            for f=1:length(elecomp)
                if isequal(aux,elecomp(f))
                    trobat=1;
                    break
                end
            end
            if (trobat==0)
                elecomp(f+1)=elecomp(f); %If the element does not exist in the elecomp, we create it
                elecomp(f)=aux; %introducing new element
            end

            % We build the matrix thanks to (f), the row index
            A(f,j)=B(f,j)+m; % The sum allows to build the matrix for species like CNN, where N=2
            B(f,j)=A(f,j);
            % We build the vector thanks to (f), the row index
            v(f)=c(f)+m*n_o(j); % The sum allows to build the matrix for species like CNN, where N=2
            c(f)=v(f);
        end
    end
    M=A;
    b=transpose(v);
end

end