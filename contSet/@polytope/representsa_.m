function [res,S] = representsa_(P,type,tol,varargin)
% representsa_ - checks if a polytope can also be represented by a
%    different set, e.g., a special case
%
% Syntax:
%    res = representsa_(P,type,tol)
%    [res,S] = representsa_(P,type,tol)
%
% Inputs:
%    P - polytope object
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

% Authors:       Mark Wetzlinger, Niklas Kochdumper, Victor Gassmann
% Written:       25-July-2023
% Last update:   01-August-2023 (MW, support fullspace comparison)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check empty object case (caution: either fullspace or empty set)
if nargout == 1
    [emptyObj,res] = representsa_emptyObject(P,type);
else
    [emptyObj,res,S] = representsa_emptyObject(P,type);
end
if emptyObj; return; end

% dimension
n = dim(P);

% init second output argument (covering all cases with res = false)
S = [];

switch type
    case 'origin'
        % quick check: is origin contained?
        res = false;
        if contains_(P,zeros(n,1),'exact',tol,0,false,false)
            % definitely not empty if the origin is contained
            P.emptySet.val = false;
            % check if only origin contained
            res = norm_(interval(P),2) <= tol;
            % set is degenerate if it's only the origin
            if res
                P.fullDim.val = false;
            end
            if nargout == 2 && res
                S = zeros(n,1);
            end
        end        

    case 'point'
        if n == 1
            V = vertices(P);
            res = size(V,2) == 1;
        else
            [fulldim,subspace] = isFullDim(P);
            res = ~fulldim && isempty(subspace);
        end
        % set is degenerate if it's only a single point
        if res
            P.fullDim.val = false;
        end
        if nargout == 2 && res
            % only one point
            % If P has a V-rep, use that one to compute the 'center'
            if P.isVRep.val
                S = center(P, 'avg');
            else
                S = center(P);
            end
        end

    case 'capsule'
        % true if 1D and bounded
        % note: also true if polytope is a bounded line
        res = n == 1 && isBounded(P);
        if nargout == 2
            throw(CORAerror('CORA:notSupported',...
                ['Conversion from polytope to ' type ' not supported.']));
        end

    case 'conHyperplane'
        if nargout == 2
            [res,S] = aux_isConHyperplane(P,tol);
        else
            res = aux_isConHyperplane(P,tol);
        end

    case 'conPolyZono'
        res = isBounded(P);
        if nargout == 2 && res
            S = conPolyZono(P);
        end

    case 'conZonotope'
        % always true for bounded polytopes
        res = isBounded(P);
        if nargout == 2 && res
            S = conZonotope(P);
        end

    case 'ellipsoid'
        % only an ellipsoid if 1D and bounded or a single point
        res = (n == 1 && isBounded(P) && ~representsa_(P,'emptySet',tol)) ...
            || representsa_(P,'point',tol);
        if nargout == 2 && res
             S = ellipsoid(P);
        end

    case 'halfspace'
        % if only a single irredundant inequality constraint given
        res = aux_isHalfspace(P,tol);

    case 'interval'
        if nargout == 2
            [res,S] = aux_isInterval(P,tol);
        else
            res = aux_isInterval(P,tol);
        end

    case 'levelSet'
        res = true;
        if nargout == 2
            S = levelSet(P);
        end

    case 'polytope'
        % obviously true
        res = true;
        if nargout == 2
            S = P;
        end

    case 'polyZonotope'
        % only true if polytope is bounded
        res = isBounded(P);
        if nargout == 2 && res
            S = polyZonotope(P);
        end

    case 'probZonotope'
        res = false;

    case 'zonoBundle'
        % zonotope bundles can represent any bounded convex set
        res = isBounded(P);
        if nargout == 2 && res
            S = zonoBundle(P);
        end

    case 'zonotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of polytope to ' type ' not supported.']));

    case 'hyperplane'
        res = isempty(P.A_.val) && size(P.Ae_.val,1) == 1;
        if res
            % hyperplanes are unbounded and non-empty
            P.bounded.val = false;
            P.emptySet.val = false;
        end

    case 'parallelotope'
        throw(CORAerror('CORA:notSupported',...
            ['Comparison of polytope to ' type ' not supported.']));

    case 'convexSet'
        res = true;

    case 'emptySet'
        res = aux_isEmptySet(P,tol);
        P.emptySet.val = res;
        % save properties, now that P is known to be the empty set
        if res
            P.bounded.val = true;
            P.fullDim.val = false;
            P.V_.val = zeros(n,0);
            P.isVRep.val = true;
            P.minVRep.val = true;
        end
        if nargout == 2 && res
            S = emptySet(n);
        end

    case 'fullspace'
        % all constraints must be trivially fulfilled, that is,
        %    A * x <= a, a >= 0, and  Ae* x = 0
        res = (P.isHRep.val ...
            && all(all(withinTol(P.A_.val,0,tol))) ...
            && all(P.b_.val > 0 | withinTol(P.b_.val,0,tol)) ...
            && all(all(withinTol(P.Ae_.val,0,tol))) ...
            && all(withinTol(P.be_.val,0,tol))) ...
            || (P.isVRep.val && (n == 1 && any(P.V_.val == -Inf) && any(P.V_.val == Inf)));
        % fullspaces are always unbounded, non-empty and full-dimensional
        if res
            P.bounded.val = false;
            P.emptySet.val = false;
            P.fullDim.val = true;
        end
        if nargout == 2 && res
            S = fullspace(n);
        end

end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_isEmptySet(P,tol)
% checks if a polytope is empty

    % check if emptiness already known
    if ~isempty(P.emptySet.val)
        res = P.emptySet.val;
    elseif P.isVRep.val
        res = size(P.V_.val,2) > 0;
    else
        res = priv_representsa_emptySet(P.A_.val,P.b_.val,...
            P.Ae_.val,P.be_.val,dim(P),tol);
    end
    P.emptySet.val = res;

end

function [res,I] = aux_isInterval(P,tol)

    % loop over all halfspaces
    % -> find those that are axis-aligned and construct an interval
    % -> gather all non-axis-aligned constraints (skip all-zero)

    % dimension and init interval
    n = dim(P);
    lb = -Inf(n,1); ub = Inf(n,1);

    all_constraints = [P.A_.val; P.Ae_.val; -P.Ae_.val];
    offsets = [P.b_.val; P.be_.val; -P.be_.val];
    nrCon = size(all_constraints,1);
    checkRedundancy = false(nrCon,1);
    
    for c=1:nrCon
        % read out constraint
        constraint = all_constraints(c,:)';

        % check if all-zero, axis-aligned, or non-axis-aligned
        if ~any(constraint)
            continue
        else
            % find non-zero index
            idx = ~withinTol(constraint,0,tol);
            if nnz(idx) == 1
                constraint_norm = vecnorm(constraint);
                % axis-aligned, save bound
                if constraint(idx) > 0
                    ub(idx) = min([offsets(c) / constraint_norm, ub(idx)]);
                else
                    lb(idx) = max([-offsets(c) / constraint_norm, lb(idx)]);
                end
            else
                % non-axis-aligned... save and check later
                checkRedundancy(c) = true;
            end
        end
    end

    % proposed interval... yet to be checked
    I = interval(lb,ub);

    % go over all non-axis-aligned constraints and check redundancy
    if any(checkRedundancy)
        for c=1:nrCon
            if checkRedundancy(c)
                % compute support function of proposed interval and compare
                % to the offset of the constraint
                val = supportFunc_(I,all_constraints(c,:)','upper');
                if val <= offsets(c)
                    % constraint is redundant
                else
                    % compute lower bound of support function to check for
                    % emptiness
                    val = supportFunc_(I,all_constraints(c,:)','lower');
                    if val > offsets(c)
                        % empty
                        res = true; I = interval.empty(n);
                        return
                    else
                        % constraint cuts through the interval... only hope
                        % now is that the polytope is empty
                        res = representsa_(P,'emptySet',tol);
                        if res
                            I = interval.empty(n);
                        end
                        return
                    end
                end
            end
        end
    end
    
    % all checks are ok
    res = true;

end

function [res,P_out] = aux_isConHyperplane(P,tol)

    % constrained hyperplane only makes sense if it is (n-1)-dimensional
    % (or empty) so that we can use it for guard intersection: this occurs
    % if and only if there is exactly one equality constraints
    if length(P.be_.val) == 1
        res = true;
        if nargout == 2
            P_out = copy(P);
        end
        return
    end

    % compute truly minimal representation and try again
    [A,b,Ae,be,empty] = priv_compact_zeros(P.A_.val,P.b_.val,P.Ae_.val,P.be_.val,tol);
    if empty
        res = false;
        P_out = [];
        return
    end
    [A,b,Ae,be] = priv_normalizeConstraints(A,b,Ae,be,'A');
    [A,b,Ae,be] = priv_compact_toEquality(A,b,Ae,be,tol);
    res = length(be) == 1;

    % instantiate 'converted' object
    if nargout == 2
        P_out = polytope(A,b,Ae,be);
        P_out = priv_copyProperties(P,P_out,'all');
    end

    % note: the set could be empty due to the inequality constraints, but
    % we skip that check for speed reasons
end

function res = aux_isHalfspace(P,tol)
% check if a polytope is only a halfspace (equivalent to a single
% inequality constraint)

%%% V representation
if P.isVRep.val
    % only 1D can be a halfspace since higher dimensions do not support
    % -Inf/+Inf vertices
    if dim(P) == 1
        % there must be either -Inf or Inf in the V representation
        % (we do not assume a minimal representation)
        signOfInfs = sign(V(isinf(V)));
        res = min(signOfInfs) == max(signOfInfs);
    else
        res = false;
    end
    return

%%% H representation
elseif P.isHRep.val

    % quick check
    if length(P.b_.val) == 1 && isempty(P.be_.val)
        res = true;
        return
    end
    
    % at least one equality constraint -> cannot be a halfspace
    if ~isempty(P.be_.val)
        res = false;
        return
    end

    % remove all-zero constraints and normalize
    [A,b,Ae,be,empty] = priv_compact_zeros(P.A_.val,P.b_.val,P.Ae_.val,P.be_.val,tol);
    
    % general check: all inequality constraints must be aligned
    [A,b,Ae,be] = priv_normalizeConstraints(A,b,Ae,be,'A');

    dotproduct = A * A';
    res = all(all(withinTol(dotproduct,1,tol)));

end

end

% ------------------------------ END OF CODE ------------------------------
