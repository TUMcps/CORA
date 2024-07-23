function val = distance(P,S)
% distance - computes the shortest distance from a polytope to another set
%    or point cloud; the cases are:
%    - one operand is the empty set: Inf
%    - the operands intersect: 0
%    - the operands do not intersect: > 0
%    - unbounded + bounded: some value < Inf
%    - unbounded + unbounded: some value < Inf  
%
% Syntax:
%    val = distance(P,S)
%
% Inputs:
%    P - polytope object
%    S - contSet object or point (cloud)
%
% Outputs:
%    val - shortest distance
%
% Example: 
%    P = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    S = polytope([ 1 0;-1 0; 0 1; 0 -1],[1;1;1;1]);
%    val = distance(P,S);
%
% Reference: MPT-Toolbox https://www.mpt3.org/
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann, Viktor Kotsev
% Written:       26-July-2021
% Last update:   07-November-2022 (added polytope-polytope case)
%                18-December-2023 (MW, add intersection check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check imput
inputArgsCheck({{P,'att','polytope','scalar'};
                {S,'att',{'ellipsoid','numeric','polytope'}}});

% set tolerance
tol = 1e-12;

% check dimensions
equalDimCheck(P,S);

% find polytope object
[P,S] = findClassArg(P,S,'polytope');

% empty set case: distance defined as infinity
if representsa(P,'emptySet') || representsa(S,'emptySet')
    val = Inf; return
end

% sets intersect (LP): shortest distance is 0
if isIntersecting_(P,S,'exact',tol)
    val = 0; return
end

% we require the halfspace representation
constraints(P);
    
% select case
if isnumeric(S)
    val = aux_distancePointCloud(P,S);

elseif isa(S,'ellipsoid')
    % call ellipsoid function
    val = distance(S,P);

elseif isa(S,'polytope')
    val = aux_distancePoly(P,S);

else
    % throw error
    throw(CORAerror('CORA:noops',P,S));

end

end


% Auxiliary functions -----------------------------------------------------

function val = aux_distancePointCloud(P,Y)

% read out dimension and number of points
[n,N] = size(Y);

% read out constraints of polytope
A = [P.A_.val; P.Ae_.val; -P.Ae_.val];
b = [P.b_.val; P.be_.val; -P.be_.val];

% init sdp problem
x = sdpvar(n,1);
C = A*x<=b;
sdpOpts = sdpsettings('verbose',0);

% multiple points or single point?
if N>1
    % construct optimizer object to increase computation speed
    y = sdpvar(n,1);
    f_obj = norm(x-y);
    dist = optimizer(C,f_obj,sdpOpts,y,f_obj);
    val = zeros(1,N);
    for i=1:N
        val(i) = dist(Y(:,i));
    end
else
    % only one point, no need for optimizer
    f_obj = norm(x-Y);
    optimize(C,f_obj,sdpOpts);
    % extract optimal point
    val = value(f_obj);
end

end

function val = aux_distancePoly(P,S)
% shortest distance between to polytopes via quadratic program

% we require the halfspace representation
constraints(S);

% put all constraints into a joint representation
problem.Aineq = blkdiag(P.A_.val, S.A_.val);
problem.bineq = [P.b_.val; S.b_.val];
problem.Aeq = blkdiag(P.Ae_.val, S.Ae_.val);
problem.beq = [P.be_.val; S.be_.val];

% read out dimension
n = dim(P);

% build the quadratic cost function (x-y)'(x-y)
nLamP = size(P.A_.val,2) - n;
nLamS = size(S.A_.val,2) - n;
HP  = diag([ones(n,1);zeros(nLamP,1)]);
HPS = blkdiag(-eye(n), zeros(nLamP, nLamS));
HS  = diag([ones(n,1);zeros(nLamS,1)]);
problem.H   = 2*[HP HPS; HPS' HS];

problem.f = zeros(size(problem.Aineq,2),1);

problem.lb = [];
problem.ub = [];

% compute using quadprog
[~,val_] = CORAquadprog(problem);
val = sqrt(val_);

end

% ------------------------------ END OF CODE ------------------------------
