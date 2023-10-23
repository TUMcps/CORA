function [res,S] = representsa_(E,type,tol,varargin)
% representsa_ - checks if a capsule can also be represented by a different
%    set, e.g., a special case
%
% Syntax:
%    res = representsa_(E,type,tol)
%    [res,S] = representsa_(E,type,tol)
%
% Inputs:
%    E - ellipsoid object
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

% Authors:       Mark Wetzlinger, Victor Gassmann
% Written:       19-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ellipsoids still allow for class arrays
if ~isscalar(E)
    if nargout == 2
        throw(CORAerror('CORA:notSupported',...
                'Second output argument not supported for ellipsoid-arrays.'));
    end

    % loop over this function
    res = false(size(E));
    for i=1:length(E)
        res(i) = representsa_(E(i),type,tol);
    end
    return
end

% check empty object case
if nargout == 1
    [empty,res] = representsa_emptyObject(E,type);
else
    [empty,res,S] = representsa_emptyObject(E,type);
end
if empty; return; end

% dimension
n = dim(E);

% init second output argument (covering all cases with res = false)
S = [];

% is the ellipsoid just a point?
isPoint = all(all(withinTol(E.Q,0,tol)));

switch type
    case 'origin'
        res = isPoint && ~isempty(E.q) && all(withinTol(E.q,0,tol));
        if nargout == 2 && res
            S = zeros(n,1);
        end

    case 'point'
        res = isPoint;
        if nargout == 2 && res
            S = E.q;
        end

    case 'capsule'
        % only if ellipsoid is 1D, a point, or a ball
        diagEQ = diag(E.Q);
        res = n == 1 || isPoint ...
            || (nnz(E.Q) == n && all(withinTol(diagEQ,diagEQ(1),tol)));
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from ellipsoid to ' type ' not supported.']));
        end

    case 'conHyperplane'
        % only a constrained hyperplane if ellipsoid is 1D or a point
        res = n == 1 || isPoint;
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from ellipsoid to ' type ' not supported.']));
        end

    case 'conPolyZono'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of ellipsoid to ' type ' not supported.']));

    case 'conZonotope'
        % only a constrained zonotope if ellipsoid is 1D or a point
        res = n == 1 || isPoint;
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from ellipsoid to ' type ' not supported.']));
        end

    case 'ellipsoid'
        % obviously true
        res = true;
        if nargout == 2
            S = E;
        end

    case 'halfspace'
        % ellipsoids cannot be unbounded
        res = false;

    case 'interval'
        % only an interval if ellipsoid is 1D or a point
        res = n == 1 || isPoint;
        if nargout == 2
            S = interval(E);
        end

    case 'levelSet'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of ellipsoid to ' type ' not supported.']));

    case 'polytope'
        % only a polytope if ellipsoid is 1D or a point
        res = n == 1 || isPoint;
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from ellipsoid to ' type ' not supported.']));
        end

    case 'polyZonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of ellipsoid to ' type ' not supported.']));

    case 'probZonotope'
        res = false;

    case 'zonoBundle'
        % only a zonotope bundle if ellipsoid is 1D or a point
        res = n == 1 || isPoint;
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from ellipsoid to ' type ' not supported.']));
        end

    case 'zonotope'
        % only a zonotope if ellipsoid is 1D or a point
        res = n == 1 || isPoint;
        if nargout == 2 && res
            S = zonotope(E);
        end

    case 'hyperplane'
        % ellipsoid cannot be unbounded (unless 1D, where hyperplane also
        % bounded)
        res = n == 1;

    case 'parallelotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of ellipsoid to ' type ' not supported.']));

    case 'emptySet'
        res = isempty(E.Q) && isempty(E.q);
        if nargout == 2 && res
            S = emptySet(n);
        end

    case 'fullspace'
        % ellipsoid cannot be unbounded
        res = false;

end

% ------------------------------ END OF CODE ------------------------------
