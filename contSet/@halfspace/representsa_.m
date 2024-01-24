function [res,S] = representsa_(hs,type,tol,varargin)
% representsa_ - checks if a halfspace can also be represented by a
%    different set, e.g., a special case
%
% Syntax:
%    res = representsa_(hs,type,tol)
%    [res,S] = representsa_(hs,type,tol)
%
% Inputs:
%    hs - halfspace object
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
    [empty,res] = representsa_emptyObject(hs,type);
else
    [empty,res,S] = representsa_emptyObject(hs,type);
end
if empty; return; end

% dimension
n = dim(hs);

% init second output argument (covering all cases with res = false)
S = [];

switch type
    case 'origin'
        res = false;

    case 'point'
        res = false;

    case 'capsule'
        res = false;

    case 'conHyperplane'
        res = false;
        
    case 'conPolyZono'
        res = false;

    case 'conZonotope'
        res = false;

    case 'ellipsoid'
        res = false;

    case 'halfspace'
        % obviously true
        res = true;
        if nargout == 2
            S = hs;
        end

    case 'interval'
        % only if normal vector is axis-aligned
        res = nnz(withinTol(hs.c,0)) >= n-1;
        if nargout == 2 && res
            S = interval(hs);
        end

    case 'levelSet'
        res = true;
        if nargout == 2 && res
            vars = sym('x',[n,1]);
            eq = hs.c' * vars - hs.d;
            S = levelSet(eq,vars,'<=');
        end

    case 'polytope'
        res = true;
        if nargout == 2
            S = polytope(hs.c,hs.d);
        end

    case 'polyZonotope'
        res = false;

    case 'probZonotope'
        res = false;

    case 'zonoBundle'
        res = false;

    case 'zonotope'
        % only empty case can also be represented by a zonotope
        res = representsa_(polytope(hs),'emptySet',tol);
        if nargout == 2 && res
            S = zonotope.empty(dim(hs));
        end

    case 'hyperplane'
        res = false;

    case 'parallelotope'
        res = false;

    case 'emptySet'
        res = representsa_(polytope(hs),'emptySet',tol);
        if nargout == 2 && res
            S = emptySet(dim(hs));
        end

    case 'fullspace'
        res = representsa_(polytope(hs),'fullspace',tol);
        if nargout == 2 && res
            S = fullspace(dim(hs));
        end

end

% ------------------------------ END OF CODE ------------------------------
