function [res,S] = representsa_(I,type,tol,varargin)
% representsa_ - checks if an interval can also be represented by a
%    different set, e.g., a special case
%
% Syntax:
%    res = representsa_(I,type,tol)
%    [res,S] = representsa_(I,type,tol)
%
% Inputs:
%    I - interval object
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
    [empty,res] = representsa_emptyObject(I,type);
else
    [empty,res,S] = representsa_emptyObject(I,type);
end
if empty; return; end

% dimension
n = dim(I);

% init second output argument (covering all cases with res = false)
S = [];

% exclude matrix interval cases
[r,c] = size(I.inf);
if ~strcmp(type,'emptySet') && r > 1 && c > 1
    throw(CORAerror('CORA:notSupported',...
        'representsa only supports vector interval objects (except type = ''emptySet'').'));
end

% is interval a point?
isPoint = all(withinTol(rad(I),0,tol));

switch type
    case 'origin'
        res = ~isempty(I.inf) ...
            && all(withinTol(I.inf,0,tol)) && all(withinTol(I.sup,0,tol));
        if nargout == 2 && res
            S = zeros(n,1);
        end

    case 'point'
        res = isPoint;
        if nargout == 2 && res
            S = center(I);
        end

    case 'capsule'
        % either 1D, a point, or at most one dimension has non-zero width
        res = n == 1 || isPoint || nnz(withinTol(rad(I),0,tol)) <= 1;
        if nargout == 2 && res
            S = capsule(I);
        end

    case 'conHyperplane'
        % only if 1D, a point, or at least one dimension has zero width
        res = n == 1 || isPoint || nnz(withinTol(rad(I),0,tol)) >= 1;
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from interval to ' type ' not supported.']));
        end
        
    case 'conPolyZono'
        res = true;
        if nargout == 2
            S = conPolyZono(I);
        end

    case 'conZonotope'
        res = true;
        if nargout == 2
            S = conZonotope(I);
        end

    case 'ellipsoid'
        % either 1D, a point, or at most one dimension has non-zero width
        res = n == 1 || isPoint || nnz(withinTol(rad(I),0,tol)) <= 1; 
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from interval to ' type ' not supported.']));
        end

    case 'halfspace'
        % only if interval is only bounded in exactly one direction
        res = nnz([reshape(isinf(I.inf),[],1); reshape(isinf(I.sup),[],1)] ) == 2*n - 1;
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from interval to ' type ' not supported.']));
        end

    case 'interval'
        % obviously true
        res = true;
        if nargout == 2
            S = I;
        end

    case 'levelSet'
        res = true;
        if nargout == 2
            S = levelSet(polytope(I));
        end

    case 'polytope'
        res = true;
        if nargout == 2
            S = polytope(I);
        end

    case 'polyZonotope'
        res = true;
        if nargout == 2
            S = polyZonotope(I);
        end

    case 'probZonotope'
        res = false;

    case 'zonoBundle'
        res = true;
        if nargout == 2
            S = zonoBundle(I);
        end

    case 'zonotope'
        res = true;
        if nargout == 2
            S = zonotope(I);
        end

    case 'hyperplane'
        % exactly one dimension has to be zero-width
        res = nnz(withinTol(rad(I),0,tol)) == 1;
        if nargout == 2 && res
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from hyperplane to ' type ' not supported.']));
        end

    case 'parallelotope'
        res = true;
        if nargout == 2
            S = zonotope(I);
        end

    case 'emptySet'
        % already handled in isemptyobject
        res = false;

    case 'fullspace'
        res = all(I.inf == -Inf) && all(I.sup == Inf);
        if nargout == 2 && res
            S = fullspace(n);
        end        

end

% ------------------------------ END OF CODE ------------------------------
