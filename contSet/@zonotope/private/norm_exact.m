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
%    Z = zonotope([1;0],[1 3 -2; 2 -1 0]);
%    norm(Z,2,'exact')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: minnorm

% Author:       Victor Gassmann
% Written:      18-September-2019
% Last update:  04-April-2023 (VG: moved YALMIP check into yalmip
%                                   "execution branch")
% Last revision:---

%------------- BEGIN CODE --------------

if ~exist('type','var')
    type = 2;
end
if type~=2
    throw(CORAerror('CORA:notSupported','Only Euclidean norm supported.'));
end

G = generators(Z);
if isempty(G)
    x = center(Z);
    val = sqrt(x'*x);
    return;
end
[~,m] = size(G);
c = center(Z);

GG = G'*G;
lmax = max(eig(GG));
M = lmax*eye(m) - GG;

persistent isMosek
if isempty(isMosek)
    isMosek = isSolverInstalled('mosek');
end

if isMosek
    % The problem to solve is 
    %%% norm(Z)^2 = max_{u\in{-1,1}^m} (c+G*u)'*(c+G*u).
    % This can be rewritten as
    %%%     -min_{u\in{-1,1}^m} u'*M*u - 2*c'*G*u - lmax*m,
    % which is the same as
    %%%     -min_{b\in{0,1}^m} 4*(b-0.5)'*M*(b-0.5) - 4*c'*G*(b-0.5)-lmax*m,
    % which is the same as
    %%%     -min_{b\in{0,1}^m} 4*||sqrtm(M)*(b-1/2*(ones(m,1)+G'*c))||_2^2
    %%%                        -c'*G*G'*c.
    % Ignoring constant terms for now, this can be rewritten as
    %%% min_{t,b} t,
    %%%     s.t.    ||sqrtm(M)*b-1/2*sqrtm(M)*(ones(m,1)+G'*c)||_2 <= t,
    %%%             0 <= b <= 1, b\in \mathbb{Z}
    
    % The second-order cone constraint can be modeled with
    %%% (t,s) \in Q^{m+1},
    % where s = sqrtm(M)*b-1/2*sqrtm(M)*(ones(m,1)+G'*c).
    % Thus, our scalar variable vector is given by x = [t;s;b]
    
    [~, res] = mosekopt('symbcon echo(0)');
    M_12 = real(sqrtm(M));
    % variables: [t;s;b] (in that order)

    % linear constraints; we need:
    %%% s = sqrtm(M)*b-1/2*sqrtm(M)*(ones(m,1)+G'*c)
    prob.a = sparse([zeros(m,1),eye(m),-M_12]);

    % lower and upper bounds are equal
    prob.blc = -1/2*M_12*(ones(m,1)+G'*c);
    prob.buc = prob.blc;

    % cones
    prob.cones.type = res.symbcon.MSK_CT_QUAD;
    prob.cones.sub = 1:m+1;
    prob.cones.subptr = 1;

    % bounds on b (0 <= b <= 1)
    prob.blx = [-inf(1,1+m),zeros(1,m)];
    prob.bux = [inf(1,1+m),ones(1,m)];

    % make b integer variable
    prob.ints.sub = 1+m+(1:m);

    % objective
    prob.c = [1,zeros(1,m+m)];

    % optimize
    [~,res] = mosekopt('minimize echo(0)',prob);

    % check if everything worked out
    if res.rcode~=0 || ~strcmp(res.sol.int.solsta,'INTEGER_OPTIMAL')
        throw(CORAerror('CORA:solverIssue','mosek'));
    end

    % extract solution
    b_sol = res.sol.int.xx(1+m+(1:m));
    u_sol = 2*(b_sol-0.5);
    x = G*u_sol+c;
    val = sqrt(-(4*(b_sol-0.5)'*M*(b_sol-0.5) - 4*c'*G*(b_sol-0.5)-lmax*m));
elseif isYalmipInstalled()
    if str2double(yalmip('version'))<20190425 % version: 25.04.2019
        throw(CORAerror('CORA:YALMIP','YALMIP version >=20190425 required'));
    end
    b = binvar(m,1);
    obj = 4*(b-0.5)'*M*(b-0.5) - 4*c'*G*(b-0.5);
    options = sdpsettings('verbose',0);
    %use gurobi for MUCH faster solve times
    %options.solver = 'bnb';
    %solve optimization problem
    optimize([], obj, options);
    val = sqrt(value(-obj + m*lmax));
    x = G*value(2*(b-0.5))+c;
    warning("YALMIP was used to model the problem - consider installing a supported solver to speed up computation...");
else
    throw(CORAerror('CORA:noSuitableSolver','integer programming'));
end

%------------- END OF CODE --------------