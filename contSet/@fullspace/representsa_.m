function [res,S] = representsa_(fs,type,tol,varargin)
% representsa_ - checks if a fullspace can also be represented by a
%    different set, e.g., a special case
%
% Syntax:
%    res = representsa_(fs,type,tol)
%    [res,S] = representsa_(fs,type,tol)
%
% Inputs:
%    fs - fullspace object
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
% Written:       24-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check empty object case
if nargout == 1
    [empty,res] = representsa_emptyObject(fs,type);
else
    [empty,res,S] = representsa_emptyObject(fs,type);
end
if empty; return; end

% dimension
n = dim(fs);

% init second output argument (covering all cases with res = false)
S = [];

switch type
    case 'origin'
        % fullspace can never be the origin
        res = false;

    case 'point'
        % fullspace can never be a single point
        res = false;

    case 'capsule'
        % capsules are bounded
        res = false;

    case 'conHyperplane'
        % constrained hyperplanes cannot cover the entire space
        res = false;

    case 'conPolyZono'
        % constrained polynomial zonotopes are bounded
        res = false;

    case 'conZonotope'
        % constrained zonotopes are bounded
        res = false;

    case 'ellipsoid'
        % ellipsoids are bounded
        res = false;

    case 'halfspace'
        % halfspace cannot cover the entire space
        res = false;

    case 'interval'
        % intervals support Inf
        res = true;
        if nargout == 2 && res
            S = interval(-Inf(n,1),Inf(n,1));
        end

    case 'levelSet'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of capsule to ' type ' not supported.']));

    case 'polytope'
        % polytopes cannot cover the entire space
        res = false;

    case 'polyZonotope'
        % polynomial zonotopes are bounded
        res = false;

    case 'probZonotope'
        res = false;

    case 'zonoBundle'
        % zonotope bundles are bounded
        res = false;

    case 'zonotope'
        % zonotopes are bounded
        res = false;

    case 'hyperplane'
        % hyperplanes cannot cover the entire space
        res = false;

    case 'parallelotope'
        % parallelotopes are bounded
        res = false;

    case 'emptySet'
        res = false;

    case 'fullspace'
        % obviously true
        res = true;
        if nargout == 2
            S = fs;
        end

end

% ------------------------------ END OF CODE ------------------------------
