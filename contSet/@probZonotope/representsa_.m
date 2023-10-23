function [res,S] = representsa_(probZ,type,tol,varargin)
% representsa_ - checks if a probabilistic zonotope can also be represented
%    by a different set, e.g., a special case
%
% Syntax:
%    res = representsa_(probZ,type,tol)
%    [res,S] = representsa_(probZ,type,tol)
%
% Inputs:
%    probZ - probZonotope object
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
    [empty,res] = representsa_emptyObject(probZ,type);
else
    [empty,res,S] = representsa_emptyObject(probZ,type);
end
if empty; return; end

% dimension
n = dim(probZ);

% init second output argument (covering all cases with res = false)
S = [];

% TODO: check if isPoint possible?

switch type
    case 'origin'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of probZonotope to ' type ' not supported.']));

    case 'point'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of probZonotope to ' type ' not supported.']));

    case 'capsule'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of probZonotope to ' type ' not supported.']));

    case 'conHyperplane'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of probZonotope to ' type ' not supported.']));
        
    case 'conPolyZono'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of probZonotope to ' type ' not supported.']));

    case 'conZonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of probZonotope to ' type ' not supported.']));

    case 'ellipsoid'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of probZonotope to ' type ' not supported.']));

    case 'halfspace'
        % probabilistic zonotope cannot be unbounded
        res = false;

    case 'interval'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of probZonotope to ' type ' not supported.']));

    case 'levelSet'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of probZonotope to ' type ' not supported.']));

    case 'polytope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of probZonotope to ' type ' not supported.']));

    case 'polyZonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of probZonotope to ' type ' not supported.']));

    case 'probZonotope'
        % obviously true
        res = true;
        if nargout == 2
            S = probZ;
        end

    case 'zonoBundle'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of probZonotope to ' type ' not supported.']));

    case 'zonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of probZonotope to ' type ' not supported.']));

    case 'hyperplane'
        % probabilistic zonotope cannot be unbounded
        res = false;

    case 'parallelotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of probZonotope to ' type ' not supported.']));

    case 'emptySet'
        res = isempty(probZ.Z) && isempty(probZ.g);

    case 'fullspace'
        % probabilistic zonotope cannot be unbounded
        res = false;

end

% ------------------------------ END OF CODE ------------------------------
