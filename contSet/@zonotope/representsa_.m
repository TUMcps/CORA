function [res,S] = representsa_(Z,type,tol,varargin)
% representsa_ - checks if a zonotope can also be represented by a
%    different set, e.g., a special case
%
% Syntax:
%    res = representsa_(Z,type,tol)
%    [res,S] = representsa_(Z,type,tol)
%
% Inputs:
%    Z - zonotope object
%    type - other set representation or 'origin', 'point', 'hyperplane'
%    tol - tolerance
%
% Outputs:
%    res - true/false
%    S - converted set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/representsa

% Authors:       Mark Wetzlinger
% Written:       19-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check empty object case
if nargout == 1
    [empty,res] = representsa_emptyObject(Z,type);
else
    [empty,res,S] = representsa_emptyObject(Z,type);
end
if empty; return; end

% dimension
n = dim(Z);

% init second output argument (covering all cases with res = false)
S = [];

% delete zero-length generators
Z = compact_(Z,'zeros',eps);
% number of generators
nrGens = size(Z.G,2);

switch type
    case 'origin'
        c = Z.c;
        res = ~isempty(c) && ( (nrGens == 0 && all(withinTol(c,0,tol))) ...
            || norm(interval(Z)) <= tol);
        if nargout == 2 && res
            S = zeros(n,1);
        end

    case 'point'
        res = nrGens == 0 || norm(interval(Z-Z.c)) <= tol;
        if nargout == 2 && res
            S = Z.c;
        end

    case 'capsule'
        % true if only one or less generators
        res = n == 1 || nrGens <= 1;
        if nargout == 2 && res
            S = capsule(Z.c,Z.G,0);
        end

    case 'conHyperplane'
        % TODO: condition?
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of zonotope to ' type ' not supported.']));
        
    case 'conPolyZono'
        res = true;
        if nargout == 2
            S = conPolyZono(Z);
        end

    case 'conZonotope'
        res = true;
        if nargout == 2
            S = conZonotope(Z);
        end

    case 'ellipsoid'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of zonotope to ' type ' not supported.']));

    case 'halfspace'
        % zonotopes cannot be unbounded
        res = false;

    case 'interval'
        res = n == 1 || aux_isInterval(Z,tol);
        if nargout == 2 && res
            S = interval(Z);
        end

    case 'levelSet'
        res = true;
        if nargout == 2
            throw(CORAerror('CORA:notSupported',...
                'Conversion from zonotope to levelSet not supported.'));
        end

    case 'polytope'
        res = true;
        if nargout == 2
            S = polytope(S);
        end

    case 'polyZonotope'
        res = true;
        if nargout == 2
            S = polyZonotope(Z);
        end

    case 'probZonotope'
        res = false;

    case 'zonoBundle'
        res = true;
        if nargout == 2
            S = zonoBundle(Z);
        end

    case 'zonotope'
        % obviously true
        res = true;
        if nargout == 2
            S = Z;
        end

    case 'hyperplane'
        % zonotopes cannot be unbounded
        res = false;

    case 'parallelotope'
        res = n == 1 || aux_isParallelotope(Z,tol);
        if nargout == 2
            S = Z;
        end

    case 'emptySet'
        % already handled in isemptyobject
        res = false;

    case 'fullspace'
        % zonotope cannot be unbounded
        res = false;

end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isInterval(Z,tol)

    res = true;
    % one-dimensional zonotopes are always intervals
    if dim(Z) == 1
        return
    end
    
    G = Z.G;
    for i=1:size(G,2)
        if nnz(~withinTol(G(:,i),0,tol)) > 1
            % two entries -> not an axis-aligned generator
            res = false;
            return
        end
    end

end

function res = aux_isParallelotope(Z,tol)

    res = true;
    
    % dimension and generators of zonotope
    n = dim(Z);
    G = Z.G;
    
    % one-dimensional zonotopes are always parallelotopes (with at least one
    % generator)
    if n == 1 && size(G,2) > 0
        return
    end
    
    % delete zero-length generators
    Z = compact_(Z,'zeros',eps);
    G = Z.G;
    
    % quick check: not enough generators
    if size(G,2) < n
        res = false; return
    end
    
    if isFullDim(Z) && size(G,2) == n
        return
    end
    
    % not a parallelotope
    res = false;

end

% ------------------------------ END OF CODE ------------------------------
