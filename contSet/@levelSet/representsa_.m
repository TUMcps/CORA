function [res,S] = representsa_(ls,type,tol,varargin)
% representsa_ - checks if a level set can also be represented by a
%    different set, e.g., a special case
%
% Syntax:
%    res = representsa_(ls,type,tol)
%    [res,S] = representsa_(ls,type,tol)
%
% Inputs:
%    ls - levelSet object
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
    [empty,res] = representsa_emptyObject(ls,type);
else
    [empty,res,S] = representsa_emptyObject(ls,type);
end
if empty; return; end

% dimension
n = dim(ls);

% init second output argument (covering all cases with res = false)
S = [];

switch type
    case 'origin'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'point'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'capsule'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'conHyperplane'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));
        
    case 'conPolyZono'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'conZonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'ellipsoid'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'halfspace'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'interval'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'levelSet'
        % obviously true
        res = true;
        if nargout == 2
            S = ls;
        end

    case 'polytope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'polyZonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'probZonotope'
        res = false;

    case 'zonoBundle'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'zonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'hyperplane'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'parallelotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of levelSet to ' type ' not supported.']));

    case 'emptySet'
        if ~iscell(ls.compOp) || length(ls.compOp) == 1
            res = false;
        else
            throw(CORAerror('CORA:notSupported',...
                'Emptiness check for multiple level sets not supported.'));
        end

    case 'fullspace'
        % level sets cannot be unbounded everywhere
        res = all(arrayfun(@(x) isempty(symvar(x)),ls.eq,'UniformOutput',true));
        if nargout == 2 && res
            S = fullspace(n);
        end

end

% ------------------------------ END OF CODE ------------------------------
