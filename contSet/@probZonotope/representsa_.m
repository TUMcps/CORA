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
    case 'probZonotope'
        % obviously true
        res = true;
        if nargout == 2
            S = probZ;
        end

    case 'hyperplane'
        % probabilistic zonotope cannot be unbounded
        res = false;

    case 'emptySet'
        res = isempty(probZ.Z) && isempty(probZ.g);

    case 'fullspace'
        % probabilistic zonotope cannot be unbounded
        res = false;

    otherwise
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of probZonotope to ' type ' not supported.']));
end

% ------------------------------ END OF CODE ------------------------------
