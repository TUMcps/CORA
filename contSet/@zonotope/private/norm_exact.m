function [val,x] = norm_exact(Z,type)
% norm_exact - computes the exact maximum norm
%
% Syntax:
%    [val,x] = norm_exact(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    val - norm value of vertex with biggest distance from the center
%    x - vertex attaining maximum norm
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: minnorm

% Authors:       Victor Gassmann
% Written:       18-September-2019
% Last update:   21-April-2023 (VG, reworked completely)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if ~isYalmipInstalled()
     throw(CORAerror('CORA:YALMIP',...
         'YALMIP must be on the MATLAB search path to use this function'));
elseif str2double(yalmip('version'))<20190425 % version: 25.04.2019
    throw(CORAerror('CORA:YALMIP','YALMIP version >=20190425 required'));
end

if ~exist('type','var')
    type = 2;
end
if type~=2
    throw(CORAerror('CORA:notSupported','Only Euclidean norm supported.'));
end

G = Z.G;
if isempty(G)
    x = Z.c;
    val = norm(x);
    return;
end
[~,m] = size(G);
c = Z.c;

GG = G'*G;
lmax = max(eig(GG));
M = lmax*eye(m) - GG;

persistent isMosek
if isempty(isMosek)
    isMosek = isSolverInstalled('mosek');
end

if isMosek
    % We want to find the vertex of the zero-centered zonotope for which
    % the distance to the center (origin) is maximized, i.e.
    %%% norm(Z)^2 = max_{x\in Z} (x-c)^T*(x-c)
    %%% norm(Z)^2 = max_{u\in{-1,1}^m} u'*G'*G*u.
    % This can be rewritten as (M = lmax*eye(m) -G'*G)
    %%%     -min_{u\in{-1,1}^m} u'*M*u - lmax*m,
    % which is the same as
    %%%     -min_{b\in{0,1}^m} 4*(b-0.5)'*M*(b-0.5) -lmax*m.
    % Ingoring the constant terms and factors, the optimizer is equal when
    % optimizing 
    %%%     min_{b\in{0,1}^m} ||sqrtm(M)*(b-0.5)||_2,
    % or equivalently
    %%% min_{t,b} t,
    %%%     s.t.    ||sqrtm(M)*b-1/2*sqrtm(M)*ones(m,1)||_2 <= t,
    %%%             0 <= b <= 1, b\in \mathbb{Z}
    
    % The second-order cone constraint can be modeled with
    %%% (t,s) \in Q^{m+1},
    % where s = sqrtm(M)*b-1/2*sqrtm(M)*ones(m,1).
    % Thus, our scalar variable vector is given by x = [t;s;b]
    
    % no linear equality constraints
    prob.a = sparse(zeros(0,1+m));
    prob.blc = [];
    prob.buc = [];

    [~, res] = mosekopt('symbcon echo(0)');
    M_sqrt = real(sqrtm(M));
    % variables: [t;b] (in that order)

    % we have the cone ||sqrtm(M)*b-1/2*sqrtm(M)*ones(m,1)||_2 <= t
    % MOSEK writes that cone as
    % (t,[sqrtm(M)*b-1/2*sqrtm(M)*ones(m,1)]_1,..) \in quadcone
    % this obviously cannot be given to the optimizer, but it can be
    % written as
    % F*[t;b] +g \in quadcone, i.e.
    prob.f = sparse([1,zeros(1,m);zeros(m,1),M_sqrt]);
    prob.g = [0;-1/2*M_sqrt*ones(m,1)];
    prob.cones = [res.symbcon.MSK_CT_QUAD, 1+m];

    % 0<=t<=inf; 0 <= b <= 1
    prob.blx = [0,zeros(1,m)];
    prob.bux = [inf,ones(1,m)];

    % make b integer variable
    prob.ints.sub = 1+(1:m);

    % objective
    prob.c = [1,zeros(1,m)];

    % optimize
    [~,res] = mosekopt('minimize echo(0)',prob);

    % check if everything worked out
    if res.rcode~=0 || ~strcmp(res.sol.int.solsta,'INTEGER_OPTIMAL')
        throw(CORAerror('CORA:solverIssue','mosek'));
    end

    % extract solution
    b_sol = res.sol.int.xx(1+(1:m));
    u_sol = 2*(b_sol-0.5);
    x = G*u_sol+c;
elseif isYalmipInstalled()
    b = binvar(m,1);
    obj = 4*(b-0.5)'*M*(b-0.5);
    options = sdpsettings('verbose',0);
    %use gurobi for MUCH faster solve times
    %options.solver = 'bnb';
    %solve optimization problem
    optimize([], obj, options);
    x = G*value(2*(b-0.5))+c;
    warning("YALMIP was used to model the problem - consider installing a supported solver to speed up computation...");
else
    throw(CORAerror('CORA:noSuitableSolver','integer programming'));
end

val = norm(x);

% ------------------------------ END OF CODE ------------------------------
