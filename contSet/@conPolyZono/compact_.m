function cPZ = compact_(cPZ,method,tol,varargin)
% compact_ - remove redundancies in the conPolyZono object representation
%
% Syntax:
%    cPZ = compact_(cPZ)
%    cPZ = compact_(cPZ,method)
%    cPZ = compact_(cPZ,method,tol)
%
% Inputs:
%    cPZ - conPolyZono object
%    method - method for redundancy removal
%             'states': remove redundancies in dependent generator matrix
%             'constraints': remove redundancies in constraint matrix
%             'exponentMatrix': remove redundancies in exponent matrix
%             'all' (default): all of the above in succession
%    tol - tolerance
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example:
%    c = [-1;0];
%    G = [1 1 0 0.5 -1 0.5; 0 0 1 1 1 0];
%    E = [1 0 0 1 2 1; 0 0 1 1 0 1; 0 0 0 1 1 1];
%    A = [0.2 -0.5 0.5 0.3];
%    b = 0.5;
%    EC = [0 1 2 0; 1 0 0 1; 0 1 0 0];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
%
%    compact(cPZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/compact, polyZonotope/compact_

% Authors:       Niklas Kochdumper
% Written:       06-November-2018
% Last update:   ---
% Last revision: 30-July-2023 (MW, rename 'compact_')

% ------------------------------ BEGIN CODE -------------------------------

    % remove all-zero independent generators (fast -> do always)
    GI = cPZ.GI;
    if ~isempty(GI)
        cPZ.GI(:,sum(abs(GI)) < tol) = [];
    end

    switch method
        case 'states'
            % remove redundancies in states
            [cPZ.c,cPZ.G,cPZ.E] = aux_removeRedundanciesGens(cPZ.c,cPZ.G,cPZ.E);

        case 'constraints'
            % remove redundancies in constraints
            if ~isempty(cPZ.A)
                [b,cPZ.A,cPZ.EC] = aux_removeRedundanciesGens(-cPZ.b,cPZ.A,cPZ.EC);
                cPZ.b = -b;
            else
                cPZ.A = []; cPZ.b = []; cPZ.EC = []; 
            end

        case 'exponentMatrix'
            % remove redudancies in the exponent matrix
            [cPZ.E,cPZ.EC,cPZ.id] = aux_removeRedundanciesExpmat(cPZ.E,cPZ.EC,cPZ.id);

        case 'all'
            % remove redundancies in states
            [cPZ.c,cPZ.G,E] = aux_removeRedundanciesGens(cPZ.c,cPZ.G,cPZ.E);
            % remove redundancies in constraints
            if ~isempty(cPZ.A)
                [b,cPZ.A,EC] = aux_removeRedundanciesGens(-cPZ.b,cPZ.A,cPZ.EC);
                cPZ.b = -b;
            else
                cPZ.A = []; cPZ.b = []; EC = []; 
            end
            % remove redudancies in the exponent matrix
            [cPZ.E,cPZ.EC,cPZ.id] = aux_removeRedundanciesExpmat(E,EC,cPZ.id);

    end

end


% Auxiliary functions -----------------------------------------------------

function [c,G,E] = aux_removeRedundanciesGens(c,G,E)
% remove redundancies in polynomial representation

    % remove redudancies in the exponent matrix
    [E,G] = removeRedundantExponents(E,G);
    
    % add all constant parts to the center
    ind = find(sum(E,1) == 0);
    
    if ~isempty(ind)
        c = c + sum(G(:,ind),2);
        G(:,ind) = [];
        E(:,ind) = [];
    end

    % catch the case where polynomial part is empty
    if isempty(G)
        c = []; G = []; E = [];
    end
end
    
function [E,EC,id] = aux_removeRedundanciesExpmat(E,EC,id)

    % remove empty rows from the exponent matrix
    E_ = [E,EC];
    temp = find(sum(E_,2) == 0);
    ind = setdiff(1:length(id),temp);
    
    E = E(ind,:); id = id(ind);
    if ~isempty(EC)
        EC = EC(ind,:); 
    end

end

% ------------------------------ END OF CODE ------------------------------
