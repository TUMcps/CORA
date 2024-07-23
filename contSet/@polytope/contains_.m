function res = contains_(P,S,type,tol,varargin)
% contains_ - determines if a polytope contains another set of points
%
% Syntax:
%    res = contains_(P,S,type,tol)
%
% Inputs:
%    P - polytope object
%    S - contSet object, numerical vector
%    tol - numerical tolerance for point in set containment
%
% Outputs:
%    res - true/false whether containment holds true
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

% Authors:       Niklas Kochdumper, Viktor Kotsev, Mark Wetzlinger
% Written:       19-November-2019
% Last update:   26-July-2021 (VG, extended to multiple points)
%                26-April-2022 (added cases for empty objects)
% Last revision: 10-July-2024 (MW, refactor)

% ------------------------------ BEGIN CODE -------------------------------

% check fullspace cases
if representsa_(P,'fullspace',0)
    % contains every set including itself
    res = true; return
elseif representsa_(S,'fullspace',0)
    % can only be contained in P if P is fullspace (would have entered the
    % if-branch above, so containment must be false)
    res = false; return
end

% point cloud in polytope containment
if isnumeric(S)
    if P.isHRep.val
        res = aux_contains_Hpoly_pointcloud(P,S,type,tol);
    else
        res = aux_contains_Vpoly_pointcloud(P,S,type,tol);
    end
    return
end

% 1D -> cheap computation of vertices (skip linear program below)
if dim(P) == 1
    try
        res = all(contains_(P,vertices(S),'exact',tol));
        return
    end
end

% nD cases
if P.isHRep.val
    % outer body in halfspace representation: check via support functions
    res = aux_contains_P_Hpoly(P,S,type,tol);
else
    % outer body in vertex representation
    res = aux_contains_P_Vpoly(P,S,type,tol);
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_contains_Hpoly_pointcloud(P,S,type,tol)
% check containment of point cloud in H-polytope, see [1, (5)]

% check inequality constraints
res_ineq = true;
if ~isempty(P.b_.val)
    offset_ineq = P.A_.val*S - P.b_.val;
    res_ineq = all(offset_ineq < tol | withinTol(offset_ineq,0,tol));
end

% check equality constraints
res_eq = true;
if ~isempty(P.be_.val)
    offset_eq = P.Ae_.val*S - P.be_.val;
    res_eq = all(offset_eq == 0 | withinTol(offset_eq,0,tol));
end

% combine checks
res = res_ineq & res_eq;

end

function res = aux_contains_Vpoly_pointcloud(P,S,type,tol)
% check containment of point cloud in V-polytope: instead of [1, (3)-(4)],
% we prefer the convhulln/Quickhull algorithm, as it is faster
% note: we output one logical value for each point in S

% read out dimension
n = dim(P);

% special method for 1D
if n == 1
    % check min/max
    V_min = min(P.V_.val);
    V_max = max(P.V_.val);
    res = S >= V_min & S <= V_max;

else
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
    
    % compute convex hull of vertices and remaining point cloud
    V_both = [V, S(:,~res)];

    try
        K = convhulln(V_both');
        % use only indices of all vertices that make up the faces of the polytope
        indices = unique(K);
        V_both_min = V_both(:,indices);
        
        % compare to original set of vertices: all vertices that are part of
        % the new convex hull are not contained in the original polytope
        for j=1:numPoints
            if ~res(j)
                res(j) = ~compareMatrices(S(:,j),V_both_min,tol,"subset");
            end
        end

    catch ME
        % likely due to degenerate point cloud... loop over points in point
        % cloud and evaluate linear program
        % min_{beta} 1
        % s.t.  v = V beta,
        %       sum_k beta_k = 1
        %       beta_k >= 0

        problem.f = zeros(numVert,1);
        problem.Aeq = [V; ones(1,numVert)];
        problem.beq = [zeros(n,1); 1];
        problem.Aineq = -eye(numVert);
        problem.bineq = zeros(numVert,1);
        problem.lb = [];
        problem.ub = [];

        for j=1:numPoints
            if ~res(j)
                % only update 'v' in LP
                problem.beq(1:n,1) = S(:,j);
                % solve linear program
                [x,fval,exitflag] = CORAlinprog(problem);
    
                % infeasible -> point not in polytope
                res(j) = exitflag ~= -2;
            end
        end
    end
end

end

function res = aux_contains_P_Hpoly(P,S,type,tol)
% containment check for any set in H-polytope

% other set is polytope in vertex representation: fast method
if isa(S,'polytope') && S.isVRep.val
    res = all(aux_contains_Hpoly_pointcloud(P,S.V_.val,type,tol));
    return
end

% generic method: check support function value along each normal vector of
% equality and inequality constraints
A = [P.A_.val; P.Ae_.val; -P.Ae_.val];
b = [P.b_.val; P.be_.val; -P.be_.val];

% additional support function options for other sets
otherOptions = {};
if isa(S,'conPolyZono') || isa(S,'polyZonotope')
    otherOptions = {'interval',8,1e-3};
end

% loop over all constraints
for i = 1:length(b)
    b_ = supportFunc_(S,A(i,:)','upper',otherOptions{:});
    if b_ > b(i) && ~withinTol(b(i),b_,tol)
        res = false;
        return 
    end
end

% all checks satisfied
res = true;

end

function res = aux_contains_P_Vpoly(P,S,type,tol)
% check containment of a set in V-polytope

% inner body is a polytope in V-representation
if isa(S,'polytope') && S.isVRep.val
    res = all(aux_contains_Vpoly_pointcloud(P,S.V_.val,type,tol));
    return
end

% if inner body is unbounded, containment has to be false since all
% V-polytopes are bounded (except for 1D, which is handled above)
if isa(S,'polytope') && ~isempty(S.bounded.val) && ~S.bounded.val
    res = false;
    return
end

% compute H-representation of outer body and check H-polytope in H-polytope
constraints(P);
res = aux_contains_P_Hpoly(P,S,type,tol);

end

% ------------------------------ END OF CODE ------------------------------
