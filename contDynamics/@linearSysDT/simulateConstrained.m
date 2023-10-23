function [t,x,ind,y] = simulateConstrained(obj,params)
% simulateConstrained - simulates a linear discrete-time system such that
% it stays within the provided reachable set; this reachable set is
% typically a backwards minmax reachable set
%
% Syntax:
%    [t,x] = simulateConstrained(obj,params)
%    [t,x,ind,y] = simulateConstrained(obj,params)
%
% Inputs:
%    obj - linearSysDT object
%    params - struct containing the parameters for the simulation
%       .tStart: initial time
%       .tFinal: final time
%       .x0: initial point
%       .U: input set
%       .W: disturbance set
%       .V: sensor noise set 
%    R - reachSet object
%
% Outputs:
%    t - time vector
%    x - state vector
%    ind - [] (argument exists only for syntax consistency)
%    y - output vector
%
% Example:
%    A = [1 2; -3 1];
%    B = [2;1];
%    dt = 1;
%    sys = linearSysDT(A,B,dt);
%
%    params.x0 = [0;0];
%    params.tFinal = 4;
%    params.u = [0.1 0.05 0.05 -0.1];
%
%    [t,x] = simulate(sys,params);
%
%    plot(x(:,1),x(:,2),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSysDT/simulate

% Authors:       Matthias Althoff
% Written:       21-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
if ~isfield(params,'tStart')
	params.tStart = 0; 
end

% set default output arguments
ind = [];
y = [];

% compute time vector and number of steps
t = (params.tStart:obj.dt:params.tFinal)';
steps = length(t)-1;

% check end of time vector
if t(end) ~= params.tFinal
    throw(CORAerror('CORA:specialError',...
        'Final time has to be a multiple of the sampling time!'));
end


% initial state
x = zeros(steps+1,obj.dim);
x(1,:) = params.x0';

% computation of output set desired / possible
comp_y = nargout == 4 && ~isempty(obj.C);

% init output
if comp_y
    y = zeros(steps+1,obj.nrOfOutputs);
end

% loop over all time steps
for i = 1:steps
    
    % sample w
    w = randPoint(params.W);
    
    % obtain current reachable set (backwards in time)
    R = params.R.timePoint.set{end-i};
    
    % compute u using linear programming
    u = aux_linProg_sol(obj, params, R, x(i,:)', w);
    
    % compute successor state
    temp = obj.A * x(i,:)' + obj.B * u + obj.c + w;
    x(i+1,:) = temp';

    % compute output
    if comp_y
        % sample v
        v = randPoint(params.V);
        % compute output
        y(i,:) = obj.C * x(i,:)' + obj.D * u + obj.k + v;
    end
    
end

% the final output is unconstrained
if comp_y
    % sample u
    u = randPoint(params.U) + params.uTrans;
    % sample v
    v = randPoint(params.V);
    % compute output
    y(i+1,:) = obj.C * x(i+1,:)' + obj.D * u + obj.k + v;
end

end


% Auxiliary functions -----------------------------------------------------

function u = aux_linProg_sol(obj, params, Z, x, w)
% obj: system object
% params: system paramters
% Z: constraint in form of a zonotope
% constraint of zonotopic reachable set: x' = c + G \beta, \beta_i \in [-1,1]
% next state: x' = Ax + Bu + w + k, u = c_u + G_u, \beta_{u,i} \in [-1,1]
% equating both states results in B G_u \beta_u - G \beta + c* = 0, 
%   c* = Ax + Bc_u + w + k - c
% objective: 
%   min \beta_u
% constraints: 
%   B G_u \beta_u - G \beta + c* = 0
%   \beta_i \in [-1,1]
%   \beta_{u,i} \in [-1,1]

% extract data
c = Z.center;
c_u = params.U.center + params.uTrans;
G = Z.generators;
G_u = params.U.generators;

% dimension and nr of generators
n_g = size(G, 2); % number of generators of Z
n_m = size(G_u, 2); % number of generators of U

% linprog solved linear programs in the form (partially using LaTex notation):
% \min_x f^T x 
% such that:
% Ax <= b \\
% A_eq x = b_eq
% x_l <= x <= x_u

% we introduce epsilon to soften the equality constraint

% define x as [\beta_u; \beta; epsilon]

% A
I = eye(n_g + n_m);
A = [I; ...
    -I];

% b
b = ones(2*(n_g + n_m),1); 
% A_eq
A_eq = [obj.B*G_u, -G];

% b_eq
b_eq = -obj.A*x - obj.B*c_u - w - obj.k + c;

% f minimizes \beta_u
f = [ones(n_m,1); zeros(n_g,1)];

% solve linear programming problem
exitflag = 0;
tol = []; % initially no tolerance
while exitflag ~= 1
    % solve linear program
    [res,~,exitflag] = linprog(f,A,b,A_eq,b_eq);
    if exitflag ~= 1
        disp('Warning! No soultion found for constrained simulation --> soften constraints')
        % enlarge tolerance
        if isempty(tol)
            tol = 1e-10;
        else
            tol = 10*tol;
        end
        % update constraints (i.e. b)
        b = [ones(n_m,1) + tol; ... % \beta_{u,i}
            ones(n_g,1); ... % \beta_i
            ones(n_m,1) + tol; ... % \beta_{u,i}
            ones(n_g,1)]; % \beta_i
    end
end

if exitflag ~= 1
    % no solution exists
    u = [];
else
    % extract \beta_u
    beta_u = res(1:n_m);
    % result
    u = c_u + G_u*beta_u;
end

end
    
% ------------------------------ END OF CODE ------------------------------
