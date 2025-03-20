function [res,cert,scaling] = contains_(P,S,method,tol,maxEval,certToggle,scalingToggle)
% contains_ - determines if a polytope contains another set of points
%
% Syntax:
%    [res,cert,scaling] = contains_(P,S,method,tol,maxEval,certToggle,scalingToggle,varargin)
%
% Inputs:
%    P - polytope object
%    S - contSet object or single point or matrix of points
%    method - method used for the containment check.
%       The available options are:
%           - 'exact': Checks for exact containment by looping over the
%               halfspaces of P.
%           - 'exact:polymax' or 'exact:venum': venum will try to check for
%               containment using the vertices of the inbody (whenever that
%               makes sense), whereas polymax tries to leverage the
%               halfspace representation of the circumbody.
%           - 'approx': If S is convex, this gives the same result as
%               'exact'. Otherwise, the convex hull of S is
%               over-approximated, and the containment check is done with
%               that approximation. This has polynomial runtime.
%               In any event, if res=1 then containment
%               is assured, whereas if res=0, S could still be contained in
%               P, but this can not be checked in polynomial time.
%    tol - tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of 
%       P will be detected as lying in P, which can be useful to 
%       counteract errors originating from floating point errors.
%    maxEval - not relevant for polytope containment
%    certToggle - if set to 'true', cert will be computed (see below).
%    scalingToggle - if set to 'true', scaling will be computed (see
%       below).
%
% Outputs:
%    res - true/false
%    cert - returns true iff the result of res could be
%           verified. For example, if res=false and cert=true, S is
%           guaranteed to not be contained in P, whereas if res=false and
%           cert=false, nothing can be deduced (S could still be
%           contained in P).
%           If res=true, then cert=true.
%           Note that computing this certification may marginally increase
%           the runtime.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(P - center(P)) + center(P) contains S.
%           Note that computing this scaling factor may significantly 
%           increase the runtime.
%           If S is convex or a point cloud, scaling is exact. If S is not
%           convex, then the returned value for scaling is an upper bound
%           of the real value.
%
% Example: 
%    P1 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = polytope([-1 -1; -1 1; 1 1;0 -1],[0;2;2;0]);
%    P3 = P2 + [2;0];
%
%    contains(P1,P2)
%    contains(P1,P3)
%
%    figure; hold on;
%    plot(P1,[1,2],'b');
%    plot(P2,[1,2],'g');
%
%    figure; hold on;
%    plot(P1,[1,2],'b');
%    plot(P3,[1,2],'r');
%
% Reference:
%    [1] M. Wetzlinger, V. Kotsev, A. Kulmburg, M. Althoff. "Implementation
%        of Polyhedral Operations in CORA 2024", ARCH'24.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains, zonotope/contains_

% Authors:       Niklas Kochdumper, Viktor Kotsev, Adrian Kulmburg, Mark Wetzlinger, Tobias Ladner
% Written:       19-November-2019
% Last update:   26-July-2021 (VG, extended to multiple points)
%                26-April-2022 (added cases for empty objects)
%                05-February-2024 (AK, added cert and scaling)
%                31-October-2024 (TL, added v-polytope/contains)
%                31-October-2024 (TL, added check if P represents a emptySet)
% Last revision: 10-July-2024 (MW, refactor)

% ------------------------------ BEGIN CODE -------------------------------

% Surprisingly, the next line is not necessarily a good idea, as it slows
% down certain instances quite a lot...
%P = compact(P);

% check fullspace cases
if representsa_(P,'fullspace',0)
    % contains every set including itself
    res = true;
    cert = true;
    scaling = 0; % This is a bit problematic here, because it's false,
    % technically speaking. The infimum is 0, but probably
    % S \not\subseteq 0 * P
    % here. There's not really a good solution here...
    return
elseif representsa_(S,'fullspace',0)
    % can only be contained in P if P is fullspace (would have entered the
    % if-branch above)
    res = false;
    cert = true;
    scaling = Inf;
    return
end

% point cloud in polytope containment
if isnumeric(S)
    % Making sure the right algorithm is used
    if ~(strcmp(method, 'exact') || strcmp(method, 'exact:venum') || ...
            strcmp(method, 'exact:polymax') || strcmp(method, 'approx'))
        throw(CORAerror('CORA:noSpecificAlg',method,P,S));
    end
    
    if P.isHRep.val || scalingToggle
        [res,cert,scaling] = aux_contains_Hpoly_pointcloud(P,S,tol,scalingToggle);
    else
        % The convex hull method can not compute the scaling
        [res,cert,scaling] = aux_contains_Vpoly_pointcloud(P,S,tol,scalingToggle);
    end

    return

end

% Now we know that S is a contSet; again, the simpler this set is, the
% better the next few algorithms might run
S = compact(S);

% First, deal with trivial cases
if representsa(S, 'emptySet')
    % Empty -> always contained
    res = true;
    cert = true;
    scaling = 0;
    return
elseif representsa(P,'emptySet')
    % S is not empty but P is -> cannot be contained
    res = false;
    cert = true;
    scaling = Inf;
    return
else
    try
        [isPoint,p] = representsa(S, 'point');
        if isPoint
            [res, cert, scaling] = contains_(P,p,'exact',tol,maxEval,certToggle,scalingToggle);
            return
        end
    catch ME
        if strcmp('CORA:notSupported', ME.identifier) || strcmp('MATLAB:maxlhs',ME.identifier)
            % If the code above returns an error either because there are
            % too many outputs, or the operation is not supported, we
            % conclude that it is not implemented, and we do nothing
        else
            % In any other case, something went wrong. Relay that
            % information.
            rethrow(ME);
        end
    end
end

% Now that the trivial cases have been dealt with, we need to choose
% which algorithm to use.

% 1D -> cheap computation of vertices (skip linear program below)
if dim(P) == 1 && ~isnumeric(S)
    [res, cert, scaling] = contains_(P,vertices(S),'exact',tol,maxEval,certToggle,scalingToggle);
    res = all(res);
    cert = true;
    scaling = max(scaling);
    return
end

% 1D -> cheap computation of vertices (skip linear program below)
if dim(P) == 1
    try
        [res,cert,scaling] = contains_(P,vertices(S),'exact',tol,maxEval,certToggle,scalingToggle);
        res = all(res);
        cert = true;
        scaling = max(scaling);
        return
    end
end

switch method
    case {'exact', 'exact:venum', 'exact:polymax'} % Exact algorithms
        [res, cert, scaling] = aux_exactParser(P, S, method, tol, maxEval, certToggle, scalingToggle);
    case 'approx' % Approximative algorithms
        % No approximative algorithm for polytopes has been implemented.
        % An alternative method is to somehow transform S into a zonotope,
        % which is approximative enough
        if ismethod(S, 'zonotope')
            [res, cert, scaling] = aux_exactParser(P, zonotope(S), 'exact', tol, maxEval, certToggle, scalingToggle);
        else
            throw(CORAerror('CORA:noSpecificAlg',method,P,S));
        end
    otherwise
        throw(CORAerror('CORA:noSpecificAlg',method,P,S));
end

end


% Auxiliary functions -----------------------------------------------------

function [res, cert, scaling] = aux_exactParser(P, S, method, tol, maxEval, certToggle, scalingToggle)
    % Chooses what algorithm to call upon to check for exact containment
    switch class(S)
        case {'conHyperplane', 'emptySet', 'fullspace', ...
                'halfspace', 'interval', 'polytope', 'conZonotope', ...
                'zonoBundle', 'zonotope', 'capsule', 'ellipsoid', ...
                'spectraShadow'}
            % All of the convex set representations
            cert = true;
            switch method
                case 'exact'
                    % Default method depends on what polytope we are
                    % dealing with
                    if P.isHRep.val
                        % outer body in halfspace representation: check via
                        % support functions
                        [res,scaling] = aux_contains_P_Hpoly(P,S,tol,scalingToggle);
                    else
                        % outer body in vertex representation
                        [res,scaling] = aux_contains_P_Vpoly(P,S,tol,scalingToggle);
                    end

                case 'exact:venum'
                    % We need to force the use of a V-polytope
                    vertices(P);
                    [res,scaling] = aux_contains_P_Vpoly(P,S,tol,scalingToggle);

                case 'exact:polymax'
                    % We need to force the use of an H-polytope
                    constraints(P);
                    [res,scaling] = aux_contains_P_Hpoly(P,S,tol,scalingToggle);
            end

        otherwise
            throw(CORAerror('CORA:noExactAlg',P,S));
    end

end


% Auxiliary functions -----------------------------------------------------

function [res, cert, scaling] = aux_contains_Hpoly_pointcloud(P,S,tol,scalingToggle)
% check containment of point cloud in H-polytope, see [1, (5)]

    % If scaling has to be computed, we need to center the polytope to the
    % origin
    if scalingToggle
        c = center(P);
        P = P - c;
        S = S - c;
    end
    % check inequality constraints
    res_ineq = true;
    if ~isempty(P.b_.val)
        offset_ineq = P.A_.val*S - P.b_.val;
        res_ineq = all(offset_ineq < tol | withinTol(offset_ineq,0,tol));
        scaling = max(P.A_.val*S./P.b_.val);
    end
    
    % check equality constraints
    res_eq = true;
    if ~isempty(P.be_.val)
        offset_eq = P.Ae_.val*S - P.be_.val;
        res_eq = all(offset_eq == 0 | withinTol(offset_eq,0,tol));
        scaling(~res_eq) = inf; % If an equality constraint is not
        % satisfied, then by convention the norm is infinite
    end
    
    % combine checks
    res = res_ineq & res_eq;
    cert = true([1 size(S, 2)]);

end

function [res,cert,scaling] = aux_contains_Vpoly_pointcloud(P,S,tol,scalingToggle)
% check containment of point cloud in V-polytope: instead of [1, (3)-(4)],
% we prefer the convhulln/Quickhull algorithm, as it is faster
% note: we output one logical value for each point in S

% Setting scaling manually, since it can not be computed using this method
scaling = zeros(1,size(S,2));
cert = true(1,size(S,2)); % The results are always correct

% read out dimension
n = dim(P);

% special method for 1D
if n == 1
    % check min/max
    V_min = min(P.V_.val);
    V_max = max(P.V_.val);
    res = S >= V_min & S <= V_max;
	return
end

% save vertices from outer body, read number of points
V = P.V;
numVert = size(V,2);
numPoints = size(S,2);

% init logical output array
res = false(1,numPoints);

% check whether points from point cloud are vertices
for j=1:numPoints
    res(j) = compareMatrices(S(:,j),V,tol,'subset');
end

% exit if all points in point cloud are vertices
if all(res)
    return
end

for i=1:numPoints
    if ~res(i)
        % add point to set of vertices at the first index
        V_added = [S(:,i) V];

        try
            K = convhulln(V_added');
            % use only indices of all vertices that make up the faces of the polytope
            indices = unique(K);
            res(i) = all(indices ~= 1);
    
        catch ME
            % likely due to degenerate point cloud... loop over points in point
            % cloud and evaluate linear program
            % min_{beta} 1
            % s.t.  v = V beta,
            %       sum_k beta_k = 1
            %       beta_k >= 0
    
            problem.f = zeros(numVert,1);
            problem.Aeq = [V; ones(1,numVert)];
            problem.beq = [S(:,i); 1];
            problem.Aineq = -eye(numVert);
            problem.bineq = zeros(numVert,1);
            problem.lb = [];
            problem.ub = [];

            % solve linear program
            [~,~,exitflag] = CORAlinprog(problem);

            % infeasible -> point not in polytope
            res(i) = exitflag ~= -2;
        end
    end

end


end

function [res,scaling] = aux_contains_P_Hpoly(P,S,tol,scalingToggle)
% containment check for any set in H-polytope

% If scaling has to be computed, first need to shift center to origin
if scalingToggle
    c = center(P);
    P = P - c;
    S = S - c;
end

% other set is polytope in vertex representation: fast method
if isa(S,'polytope') && S.isVRep.val
    [res_pc,~,scaling_pc] = aux_contains_Hpoly_pointcloud(P,S.V_.val,tol,scalingToggle);
    res = all(res_pc);
    scaling = max(scaling_pc);
    return
end

% generic method: check support function value along each normal vector of
% equality and inequality constraints
[A,b] = priv_equalityToInequality(P.A_.val,P.b_.val,P.Ae_.val,P.be_.val);

% additional support function options for other sets
otherOptions = {};
if isa(S,'conPolyZono') || isa(S,'polyZonotope')
    otherOptions = {'interval',8,1e-3};
end

% loop over all constraints
scaling = 0;
res = true;

for i = 1:length(b)
    b_ = supportFunc_(S,A(i,:)','upper',otherOptions{:});
    if b_ > b(i) && ~withinTol(b(i),b_,tol)
        res = false;
        if ~scalingToggle
            % If scaling does not have to be computed, we can stop here, to
            % speed up the procedure
            return
        else
            % Otherwise, we need to extract the scaling.
            scaling_part = b_/b(i);
            % A couple of things can happen here because of degeneracy:
            % If scaling_part = inf this means the inbody is unbounded, but
            % all is good and scaling has been correctly computed.
            % On the other hand, if scaling_part is NaN, we can set it to
            % zero, as it means the inbody and circumbody are both
            % degenerate in that direction.
            if isnan(scaling_part)
                scaling_part = 0;
            end
            scaling = max([scaling scaling_part]);
        end
    end
end

end

function [res, scaling] = aux_contains_P_Vpoly(P,S,tol,scalingToggle)
% check containment of a set in V-polytope

% inner body is a polytope in V-representation
if isa(S,'polytope') && S.isVRep.val
    if scalingToggle
        constraints(P); % Need to compute Hrep, otherwise we can't compute
                        % the scaling
        [res_pc,~,scaling_pc] = aux_contains_Hpoly_pointcloud(P,S.V_.val,tol,scalingToggle);
    else
        [res_pc,~,scaling_pc] = aux_contains_Vpoly_pointcloud(P,S.V_.val,tol,scalingToggle);
    end
        res = all(res_pc);
        scaling = max(scaling_pc);
        return
    end

% if inner body is unbounded, containment has to be false since all
% V-polytopes are bounded (except for 1D, which is handled above)
if isa(S,'polytope') && ~isempty(S.bounded.val) && ~S.bounded.val
    res = false;
    scaling = inf;
    return
end

% compute H-representation of outer body and check H-polytope in H-polytope
constraints(P);
[res, scaling] = aux_contains_P_Hpoly(P,S,tol,scalingToggle);
end

% ------------------------------ END OF CODE ------------------------------
