function [res,S] = representsa_(hyp,type,tol,varargin)
% representsa_ - checks if a constrained hyperplane can also be represented
%    by a different set, e.g., a special case
%
% Syntax:
%    res = representsa_(hyp,type,tol)
%    [res,S] = representsa_(hyp,type,tol)
%
% Inputs:
%    hyp - conHyperplane object
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

% check empty object case
if nargout == 1
    [empty,res] = representsa_emptyObject(hyp,type);
else
    [empty,res,S] = representsa_emptyObject(hyp,type);
end
if empty; return; end

% dimension
n = dim(hyp);

% init second output argument (covering all cases with res = false)
S = [];

switch type
    case 'origin'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conHyperplane to ' type ' not supported.']));

    case 'point'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conHyperplane to ' type ' not supported.']));

    case 'capsule'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conHyperplane to ' type ' not supported.']));

    case 'conHyperplane'
        % obviously true
        res = true;
        if nargout == 2
            S = hyp;
        end

    case 'conPolyZono'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conHyperplane to ' type ' not supported.']));

    case 'conZonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conHyperplane to ' type ' not supported.']));

    case 'ellipsoid'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conHyperplane to ' type ' not supported.']));

    case 'halfspace'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conHyperplane to ' type ' not supported.']));

    case 'interval'
        % only if normal vector of halfspace axis-aligned
        res = n == 1 || nnz(withinTol(hyp.a',0,tol)) >= n-1;
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from conHyperplane to ' type ' not supported.']));
        end

    case 'levelSet'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conHyperplane to ' type ' not supported.']));

    case 'polytope'
        res = true;
        if nargout == 2
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from conHyperplane to ' type ' not supported.']));
        end

    case 'polyZonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conHyperplane to ' type ' not supported.']));

    case 'probZonotope'
        res = false;

    case 'zonoBundle'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conHyperplane to ' type ' not supported.']));

    case 'zonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conHyperplane to ' type ' not supported.']));

    case 'hyperplane'
        res = n == 1 || aux_isHyperplane(hyp,tol);
        if nargout == 2 && res
            S = hyp;
        end

    case 'parallelotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of conHyperplane to ' type ' not supported.']));

    case 'emptySet'
        res = representsa_(polytope(hyp.C,hyp.d,hyp.a,hyp.b),'emptySet',tol);
        if nargout == 2 && res
            S = emptySet(n);
        end

    case 'fullspace'
        res = representsa_(polytope(hyp.C,hyp.d,hyp.a,hyp.b),'fullspace',tol);
        if nargout == 2 && res
            S = fullspace(n);
        end


end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isHyperplane(hyp,tol)

    % always true if conHyperplane is 1D and consistent (checked in constructor)
    if dim(hyp)==1
        res = true;
        return;
    end
    if isempty(hyp.C) || (all(all(hyp.C==0)) && all(hyp.d>=0))
        res = true;
        return;
    end
    
    % check if C is unbounded for all x on the hyperplane by checking if
    % support function is unbounded for all 2*(n-1) directions ('upper' and 'lower')
    c = hyp.a'/norm(hyp.a');
    n = dim(hyp);
    % null space has exactly n-1 vectors
    B = null(c');
    for i=1:n-1
        if supportFunc_(hyp,B(:,1),'upper') < Inf || ...
           supportFunc_(hyp,B(:,1),'lower') > -Inf
            res = false;
            return;
        end
    end
    
    res = true;

end

% ------------------------------ END OF CODE ------------------------------
