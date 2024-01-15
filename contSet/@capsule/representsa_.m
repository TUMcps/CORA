function [res,S] = representsa_(C,type,tol,varargin)
% representsa_ - checks if a capsule can also be represented by a different
%    set, e.g., a special case
%
% Syntax:
%    res = representsa_(C,type,tol)
%    [res,S] = representsa_(C,type,tol)
%
% Inputs:
%    C - capsule object
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
    [empty,res] = representsa_emptyObject(C,type);
else
    [empty,res,S] = representsa_emptyObject(C,type);
end
if empty; return; end

% dimension
n = dim(C);

% init second output argument (covering all cases with res = false)
S = [];

% is the capsule just a point?
isPoint = ~isempty(C.c) && all(withinTol(C.g,zeros(n,1),tol)) ...
    && withinTol(C.r,0,tol);

switch type
    case 'origin'
        res = isPoint && all(withinTol(C.c,0,tol));
        if nargout == 2 && res
            S = zeros(n,1);
        end

    case 'point'
        res = isPoint;
        if nargout == 2 && res
            S = C.c;
        end

    case 'capsule'
        % obviously true
        res = true;
        if nargout == 2
            S = C;
        end

    case 'conHyperplane'
        % only a constrained hyperplane if 1D or no radius
        res = n == 1 || withinTol(C.r,0,tol);
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from capsule to ' type ' not supported.']));
        end

    case 'conPolyZono'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of capsule to ' type ' not supported.']));

    case 'conZonotope'
        % only a constrained zonotope if 1D or no radius
        res = n == 1 || withinTol(C.r,0,tol);
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from capsule to ' type ' not supported.']));
        end

    case 'ellipsoid'
        % only an ellipsoid if 1D or no generator or generator, but no
        % radius
        res = n == 1 || all(withinTol(C.g,0,tol)) || withinTol(C.r,0,tol);
        if nargout == 2 && res
             S = ellipsoid(C);
        end

    case 'halfspace'
        % capsules cannot be unbounded
        res = false;

    case 'interval'
        % only an interval if no radius and generator is axis-aligned
        res = n == 1 || ...
            (withinTol(C.r,0,tol) && nnz(withinTol(C.g,0,tol)) >= n-1);
        if nargout == 2 && res
            S = interval(C);
        end

    case 'levelSet'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of capsule to ' type ' not supported.']));

    case 'polytope'
        % only a polytope if 1D or no radius
        res = n == 1 || withinTol(C.r,0,tol);
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from capsule to ' type ' not supported.']));
        end

    case 'polyZonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of capsule to ' type ' not supported.']));

    case 'probZonotope'
        res = false;

    case 'zonoBundle'
        % only a zonotope bundle if 1D or no radius
        res = n == 1 || withinTol(C.r,0,tol);
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from capsule to ' type ' not supported.']));
        end

    case 'zonotope'
        % only a zonotope if 1D or no radius
        res = n == 1 || withinTol(C.r,0,tol);
        if nargout == 2 && res
            S = zonotope(C);
        end

    case 'hyperplane'
        % capsule cannot be unbounded (except 1D where hyperplane also
        % bounded)
        res = n == 1;

    case 'parallelotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of capsule to ' type ' not supported.']));

    case 'emptySet'
        res = isempty(center(C)) || isempty(C.g) || isempty(C.r);
        if nargout == 2 && res
            S = emptySet(n);
        end

    case 'fullspace'
        % capsule cannot be unbounded
        res = false;

end

% ------------------------------ END OF CODE ------------------------------
