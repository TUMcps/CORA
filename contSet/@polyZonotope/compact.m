function res = compact(pZ)
% compact - remove redundancies in the representation of a poly. zonotope 
%
% Syntax:  
%    res = compact(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    res - polyZonotope object without redundancies
%
% Example: 
%    pZ = polyZonotope([1;2],[1 3 1 -1 0;0 1 1 1 2], ...
%                      [],[1 0 1 0 1;0 1 2 0 2])
%    compact(pZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conPolyZono/compact, zonotope/deleteZeros

% Author:       Niklas Kochdumper
% Written:      20-January-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
    % remove redudancies in the exponent matrix
    [expMat,G] = removeRedundantExponents(pZ.expMat,pZ.G);
    
    % add all constant parts to the center
    ind = find(sum(expMat,1) == 0);
    
    if ~isempty(ind)
        c = pZ.c + sum(G(:,ind),2);
        G(:,ind) = [];
        expMat(:,ind) = [];
    else
        c = pZ.c;
    end
    
    % remove empty generators
    Grest = pZ.Grest;
    if ~isempty(Grest)
        Grest(:,sum(abs(Grest)) < eps) = [];
    end
    
    % remove empty rows from the exponent matrix
    id = pZ.id;
    ind = find(sum(expMat,2) == 0);
    if ~isempty(ind)
        expMat(ind,:) = [];
        id(ind) = [];
    end
    
    % construct compacted polynomial zonotope
    res = polyZonotope(c,G,Grest,expMat,id);
end
   
%------------- END OF CODE --------------