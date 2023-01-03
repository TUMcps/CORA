function res = isemptyobject(S)
% isemptyobject - checks whether a contSet object contains any information
%    at all; consequently, the set is equivalent to the empty set 
%
% Syntax:  
%    res = isemptyobject(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    res - true/false
%
% Example: 
%    E = ellipsoid([1,0;0,1]);
%    isemptyobject(E);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: isempty

% Author:       Mark Wetzlinger
% Written:      03-June-2022
% Last update:  18-August-2022 (MW, extend to class arrays)
% Last revision:---

%------------- BEGIN CODE --------------

% read class of contSet object
classname = class(S);

% differentiate between interval, taylm, affine, zoo (range-bounding) and
% more 'standard' contSet classes
if any(strcmp(classname,{'affine','interval','taylm','zoo'}))
    % sets displayed as nx1 objects are treated as 1x1 objects
    res = checkIfEmpty(S);
else % all other contSet classes (can be nx1 arrays)
    res = false(size(S));
    % loop over class-arrays
    for i=1:size(S,1)
        for j=1:size(S,2)
                res(i,j) = checkIfEmpty(S(i,j));
        end
    end
end

end

% Auxiliary function ------------------------------------------------------

function res = checkIfEmpty(S)
% properties must be same as in fully-empty object of same class
% this is a rather brute-force version that needs to be updated
% semi-regularly, but for efficiency reasons there is no real other option

classname = class(S);

% different checks depending on class
switch classname
    case 'affine'
        % no empty object allowed
        res = false;

    case 'capsule'
        res = isnumeric(S.c) && isempty(S.c) ...
            && isnumeric(S.g) && isempty(S.g) ...
            && isnumeric(S.r) && isscalar(S.r) && S.r == 0;

    case 'conHyperplane'
        res = isa(S.h,'halfspace') ...
            && isnumeric(S.h.c) && isempty(S.h.c) ...
            && isnumeric(S.h.d) && isscalar(S.h.d) && S.h.d == 0 ...
            && isnumeric(S.C) && isempty(S.C) ...
            && isnumeric(S.d) && isscalar(S.d) && S.d == 0;

    case 'conPolyZono'
        res = isnumeric(S.c) && isempty(S.c) ...
            && isnumeric(S.G) && isempty(S.G) ...
            && isnumeric(S.expMat) && isempty(S.expMat) ...
            && isnumeric(S.A) && isempty(S.A) ...
            && isnumeric(S.b) && isempty(S.b) ...
            && isnumeric(S.expMat_) && isempty(S.expMat_) ...
            && isnumeric(S.Grest) && isempty(S.Grest) ...
            && isnumeric(S.id) && isempty(S.id);

    case 'conZonotope'
        res = isnumeric(S.Z) && isempty(S.Z) ...
            && isnumeric(S.A) && isempty(S.A) ...
            && isnumeric(S.b) && isempty(S.b) ...
            && isnumeric(S.ksi) && isempty(S.ksi) ...
            && isnumeric(S.R) && isempty(S.R);

    case 'ellipsoid'
        res = isnumeric(S.Q) && isempty(S.Q) ...
            && isnumeric(S.q) && isempty(S.q) ...
            && isnumeric(S.TOL) && isscalar(S.TOL) && S.TOL == 1e-6;

    case 'halfspace'
        res = isnumeric(S.c) && isempty(S.c) ...
            && isnumeric(S.d) && isscalar(S.d) && S.d == 0;

    case 'interval'
        res = isnumeric(S.inf) && isempty(S.inf) ...
            && isnumeric(S.sup) && isempty(S.sup);

    case 'levelSet'
        res = isnumeric(S.eq) && isempty(S.eq) ...
            && isnumeric(S.vars) && isempty(S.vars) ...
            && isnumeric(S.compOp) && isempty(S.compOp) ...
            && isnumeric(S.funHan) && isempty(S.funHan) ...
            && isnumeric(S.der) && isempty(S.der) ...
            && isnumeric(S.dim) && isempty(S.dim) ...
            && isnumeric(S.solved) && isempty(S.solved) ...
            && islogical(S.solvable) && ~S.solvable;

    case 'mptPolytope'
        % wrapper for Polyhedron class
        res = isnumeric(S.P) && isempty(S.P) ...
            && isnumeric(S.halfspace) && isempty(S.halfspace);

    case 'Polyhedron'
        % from mpt toolbox
        res = islogical(S.irredundantVRep) && ~S.irredundantVRep ...
            && islogical(S.irredundantHRep) && ~S.irredundantHRep ...
            && islogical(S.hasHRep) && ~S.hasHRep ...
            && islogical(S.hasVRep) && ~S.hasVRep ...
            && isnumeric(S.A) && isempty(S.A) ...
            && isnumeric(S.b) && isempty(S.b) ...
            && isnumeric(S.Ae) && isempty(S.Ae) ...
            && isnumeric(S.be) && isempty(S.be) ...
            && isnumeric(S.H) && isempty(S.H) ...
            && isnumeric(S.He) && isempty(S.He) ...
            && isnumeric(S.R) && isempty(S.R) ...
            && isnumeric(S.V) && isempty(S.V) ...
            && isnumeric(S.Dim) && isscalar(S.Dim) && S.Dim == 0 ...
            && isnumeric(S.Data) && isempty(S.Data);

    case 'polyZonotope'
        res = isnumeric(S.c) && isempty(S.c) ...
            && isnumeric(S.G) && isempty(S.G) ...
            && isnumeric(S.Grest) && isempty(S.Grest) ...
            && isnumeric(S.expMat) && isempty(S.expMat) ...
            && isnumeric(S.id) && isempty(S.id);

    case 'probZonotope'
        res = isnumeric(S.Z) && isempty(S.Z) ...
            && isnumeric(S.g) && isempty(S.g) ...
            && isnumeric(S.cov) && isempty(S.cov) ...
            && islogical(S.gauss) && ~S.gauss ...
            && isnumeric(S.gamma) && isscalar(S.gamma) && S.gamma == 2;

    case 'taylm'
        res = isnumeric(S.coefficients) && isscalar(S.coefficients) && S.coefficients == 0 ...
            && isnumeric(S.monomials) && isempty(S.monomials) ...
            && isa(S.remainder,'interval') ...
            && isscalar(S.remainder.inf) && isscalar(S.remainder.sup) ...
            && S.remainder.inf == 1 && S.remainder.sup == 1 ...
            && iscell(S.names_of_var) && isempty(S.names_of_var) ...
            && isnumeric(S.max_order) && isscalar(S.max_order) && S.max_order == 6 ...
            && ischar(S.opt_method) && strcmp(S.opt_method,'int') ...
            && isnumeric(S.eps) && isscalar(S.eps) && S.eps == 1e-3 ...
            && isnumeric(S.tolerance) && isscalar(S.tolerance) && S.tolerance == 1e-8;

    case 'zonoBundle'
        res = iscell(S.Z) && isempty(S.Z) ...
            && isnumeric(S.parallelSets) && isscalar(S.parallelSets) && S.parallelSets == 0;

    case 'zonotope'
        res = isnumeric(S.Z) && isempty(S.Z) ...
            && isnumeric(S.halfspace) && isempty(S.halfspace);

    case 'zoo'
        res = false;

end

end

function res = checkIfEmpty_old(S,S_empty,props)
% old method based on instantiated empty object (which takes a long time,
% therefore replaced by more direct albeit less sustainable solution)
% this function is kept for debugging purposes;
% note: compute other two input arguments by
%    S_empty = eval([classname '();']);
%    props = properties(S);

% properties must be same as in fully-empty object of same class

    % init return value
    res = true;

    nrProps = length(props);
    for p=1:nrProps
    
        % property of contSet object
        prop = S.(props{p});
    
        % corresponding property of fully-empty set
        prop_empty = S_empty.(props{p});
    
        % type of property
        type = class(prop);
    
        % check depends on type
        switch type
    
            case 'logical'
                if prop ~= prop_empty
                    res = false; return
                end
    
            case {'double','cell','sym'}
                if ~isequal(prop,prop_empty)
                    res = false; return
                end
    
            case {'char','string'}
                if ~strcmp(prop,prop_empty)
                    res = false; return
                end
    
            case 'Polyhedron'
                % for empty mptPolytope objects, we have [] instead of a
                % Polyhedron object, thus it's always false
                res = false; return
                % caution: this check will be changed once the mptPolytope
                % class is substituted
    
            otherwise
                % check if the property is a contSet object (e.g., halfspace
                % object in conHyperplane object)
                if isa(prop,'contSet')
                    % recursive call
                    if ~isemptyobject(prop)
                        res = false; return
                    end
                else
                    throw(CORAerror('CORA:specialError','Unknown data type'));
                end
        end
        
    end

end

%------------- END OF CODE --------------