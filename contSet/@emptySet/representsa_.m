function [res,S] = representsa_(O,type,tol,varargin)
% representsa_ - checks if an empty set can also be represented by a
%    different set, e.g., a special case
%
% Syntax:
%    res = representsa_(O,type,tol)
%    [res,S] = representsa_(O,type,tol)
%
% Inputs:
%    O - emptySet object
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
% Written:       25-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check empty object case
if nargout == 1
    [empty,res] = representsa_emptyObject(O,type);
else
    [empty,res,S] = representsa_emptyObject(O,type);
end
if empty; return; end

% dimension
n = dim(O);

% init second output argument (covering all cases with res = false)
S = [];

switch type
    case 'origin'
        % empty set can never be the origin
        res = false;

    case 'point'
        % empty set can never be a point
        res = false;

        % all contSet can be empty
    case 'capsule'
        res = true;
        if nargout == 2 % call static function
            S = capsule.empty(dim(O));
        end

    case 'conHyperplane'
        res = true;
        if nargout == 2 % call static function
            S = polytope.empty(dim(O));
        end

    case 'conPolyZono'
        res = true;
        if nargout == 2 % call static function
            S = conPolyZono.empty(dim(O));
        end

    case 'conZonotope'
        res = true;
        if nargout == 2 % call static function
            S = conZonotope.empty(dim(O));
        end

    case 'ellipsoid'
        res = true;
        if nargout == 2 % call static function
            S = ellipsoid.empty(dim(O));
        end

    case 'halfspace'
        % empty set is always empty
        res = false;

    case 'interval'
        res = true;
        if nargout == 2 % call static function
            S = interval.empty(dim(O));
        end

    case 'levelSet' % not supported
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of capsule to ' type ' not supported.']));

    case 'polytope'
        res = true;
        if nargout == 2 % call static function
            S = polytope.empty(dim(O));
        end

    case 'polyZonotope'
        res = true;
        if nargout == 2 % call static function
            S = polyZonotope.empty(dim(O));
        end

    case 'probZonotope'
        % cannot be true
        res = false;

    case 'zonoBundle'
        res = true;
        if nargout == 2 % call static function
            S = zonoBundle.empty(dim(O));
        end

    case 'zonotope'
        res = true;
        if nargout == 2
            S = zonotope.empty(dim(O));
        end

    case 'hyperplane'
        % hyperplanes cannot be empty
        res = false;

    case 'parallelotope'
        res = false;

    case 'convexSet'
        res = true;

    case 'emptySet'
        % obviously true
        res = true;
        if nargout == 2
            S = O;
        end

    case 'fullspace'
        res = false;

end

% ------------------------------ END OF CODE ------------------------------
