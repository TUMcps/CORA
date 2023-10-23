function pZ = compact_(pZ,method,tol,varargin)
% compact_ - removes redundancies in the representation of a polynomial
%    zonotope, such as zero-length generators or unused dependent factors
%
% Syntax:
%    pZ = compact_(pZ)
%    pZ = compact_(pZ,method)
%    pZ = compact_(pZ,method,tol)
%
% Inputs:
%    pZ - polyZonotope object
%    method - method for redundancy removal
%             'states': remove redundancies in dependent generator matrix
%             'exponentMatrix': remove redundancies in exponent matrix
%             'all' (default): all of the above in succession
%    tol - tolerance
%
% Outputs:
%    pZ - polyZonotope object without redundancies
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
% See also: contSet/compact, conPolyZono/compact_, zonotope/compact_

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       20-January-2020
% Last update:   ---
% Last revision: 30-July-2023 (MW, restructure, merge with deleteZeros)

% ------------------------------ BEGIN CODE -------------------------------

switch method
    case 'states'
        % remove empty generators (fast -> always)
        [pZ.G,pZ.GI,pZ.E,pZ.id] = ...
            aux_removeRedundanciesGI(pZ.G,pZ.GI,pZ.E,pZ.id,tol);

    case 'exponentMatrix'
        % remove redudancies in the exponent matrix
        [pZ.c,pZ.G,pZ.E,pZ.id] = ...
            aux_removeRedundanciesE(pZ.c,pZ.G,pZ.E,pZ.id);

    case 'all'
        % remove empty generators (fast -> always)
        [pZ.G,pZ.GI,pZ.E,pZ.id] = ...
            aux_removeRedundanciesGI(pZ.G,pZ.GI,pZ.E,pZ.id,tol);
        % remove redudancies in the exponent matrix
        [pZ.c,pZ.G,pZ.E,pZ.id] = ...
            aux_removeRedundanciesE(pZ.c,pZ.G,pZ.E,pZ.id);

end

end


% Auxiliary functions -----------------------------------------------------

function [G,GI,E,id] = aux_removeRedundanciesGI(G,GI,E,id,tol)
    
    % indices with non-zero generators
    idxD = any(G,1);
    idxI = any(GI,1);
    
    % if all non-zero, skip
    if ~( all(idxD) && all(idxI) ) 
        % delete zero generators
        G = G(:,idxD);
        E = E(:,idxD);
        GI = GI(:,idxI);
        
        % delete zero exponents
        idxE = any(E,2);
        if ~all(idxE)
            E = E(idxE,:);
            id = id(idxE);
        end
    else
        return
    end

end

function [c,G,E,id] = aux_removeRedundanciesE(c,G,E,id)

    [E,G] = removeRedundantExponents(E,G);
        
    % add all constant parts to the center
    ind = find(sum(E,1) == 0);
    
    if ~isempty(ind)
        c = c + sum(G(:,ind),2);
        G(:,ind) = [];
        E(:,ind) = [];
    end
    
    % remove empty rows from the exponent matrix
    ind = find(sum(E,2) == 0);
    if ~isempty(ind)
        E(ind,:) = [];
        id(ind) = [];
    end

end

% ------------------------------ END OF CODE ------------------------------
