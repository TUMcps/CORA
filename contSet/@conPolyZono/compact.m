function cPZ = compact(cPZ)
% compact - remove redundancies in the conPolyZono representation
%
% Syntax:  
%    cPZ = compact(cPZ)
%
% Inputs:
%    cPZ - conPolyZono object
%
% Outputs:
%    cPZ - corresponding conPolyZono object without redundancies
%
% Example:
%    c = [-1;0];
%    G = [1 1 0 0.5 -1 0.5; 0 0 1 1 1 0];
%    expMat = [1 0 0 1 2 1; 0 0 1 1 0 1; 0 0 0 1 1 1];
%    A = [0.2 -0.5 0.5 0.3];
%    b = 0.5;
%    expMat_ = [0 1 2 0; 1 0 0 1; 0 1 0 0];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%
%    compact(cPZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/compact

% Author:       Niklas Kochdumper
% Written:      06-November-2018
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE -------------

    % remove redundancies in states
    [c,G,expMat] = removeRedundancies(cPZ.c,cPZ.G,cPZ.expMat);
    
    % remove redundancies in constraints
    if ~isempty(cPZ.A)
        [b,A,expMat_] = removeRedundancies(-cPZ.b,cPZ.A,cPZ.expMat_);
        b = -b;
    else
        A = []; b = []; expMat_ = []; 
    end
    
    % remove empty generators
    Grest = cPZ.Grest;
    if ~isempty(Grest)
        Grest(:,sum(abs(Grest)) < eps) = [];
    end

    % remove empty rows from the exponent matrix
    E = [expMat,expMat_];
    temp = find(sum(E,2) == 0);
    ind = setdiff(1:length(cPZ.id),temp);
    
    expMat = expMat(ind,:); id = cPZ.id(ind);
    if ~isempty(expMat_)
        expMat_ = expMat_(ind,:); 
    end
    
    % construct compacted constrained polynomial zonotope
    cPZ = conPolyZono(c,G,expMat,A,b,expMat_,Grest,id);    
end


% Auxiliary Functions -----------------------------------------------------

function [c,G,expMat] = removeRedundancies(c,G,expMat)
% remove redundancies in polynomial representation

    % remove redudancies in the exponent matrix
    [expMat,G] = removeRedundantExponents(expMat,G);
    
    % add all constant parts to the center
    ind = find(sum(expMat,1) == 0);
    
    if ~isempty(ind)
        c = c + sum(G(:,ind),2);
        G(:,ind) = [];
        expMat(:,ind) = [];
    end

    % catch the case where polynomial part is empty
    if isempty(G)
        c = []; G = []; expMat = [];
    end
end
    
%------------- END OF CODE --------------