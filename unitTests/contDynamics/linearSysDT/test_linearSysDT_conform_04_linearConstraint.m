function res = test_linearSysDT_conform_04_linearConstraint
% test_linearSysDT_conform_04_linearConstraint - unit test  
%   for checking the computation of the linear constraints for conformance 
%   synthesis of linear discrete-time systems according to 
%   [1] and [2].
%
% Syntax:
%    res = test_linearSysDT_conform_04_linearConstraint
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] M. Althoff. Checking and Establishing Reachset Conformance in
%        CORA 2023. In Proc. of the 10th International Workshop on 
%        Applied Verification of Continuous and Hybrid Systems, 2023.
%    [2] Liu et al., "Guarantees for Real Robotic Systems: Unifying Formal
%        Controller Synthesis and Reachset-Conformant Identification", 202x.

% Authors:       Matthias Althoff
% Written:       12-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init partial results
resPartial = [];

%% System dynamics
% create pedestrian model (see Lecture "Formal Methods for Cyber-Physical
% Systems - Conformance Checking")
% sample time
dt = 0.4;
% discrete time system matrix from continuous system matrix Ac
Ac = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
A = expm(Ac*dt);
B = [];
C = [1 0 0 0; 0 1 0 0];
D = [];

% instantiate linear discrete-time dynamics
pedestrian = linearSysDT(A,B,[],C,D,dt);
params.tFinal = 4*dt;

% the system has known uncertainties
params.W = zonotope(rand(4,5)-0.5); % disturbance 
params.V = zonotope(rand(2,3)-0.5); % measurement uncertainty 
params.R0 = zonotope(rand(4,3)-0.5); % initial state uncertainty; initial state uncertainty has no effect for linear systems

%% Compute reachable set
options.zonotopeOrder = inf;
R = reach(pedestrian, params, options);

%% create test cases as vertices of reachable set
for iStep = 1:length(R.timePoint.time)
    union_y_a{iStep} = vertices(R.timePoint.set{iStep}).';
end

% change final time for analysis
params.tFinal = 3*dt;

%% compute linear constraint A*z <= b
[lin_A, lin_b, G_total_saved, Gamma_saved, N_saved, G_tilde_saved] = aux_linearConstraints(pedestrian,union_y_a, params);

%% assemble center vector
c = [center(params.R0); center(params.W); center(params.V)];

%% compare reachable sets assembled from G_total_saved, Gamma_saved
for iStep = 1:length(R.timePoint.time)
    % assemble reachable set from computing linear constraints
    R_assembled = zonotope([Gamma_saved{iStep}*c,G_total_saved{iStep}]);
    % enlarged R_assembled should contain R
    R_enlarged = enlarge(R_assembled,1.001);
    resPartial(end+1) = contains(R_enlarged, R.timePoint.set{iStep});
    % enlarged R should contain R_assembled
    R_enlarged = enlarge(R.timePoint.set{iStep},1.001);
    resPartial(end+1) = contains(R_enlarged, R_assembled);
end

%% compare reachable sets assembled from halfspace representation
for iStep = 1:length(R.timePoint.time)
    % assemble reachable set from computing linear constraints
    N = N_saved{iStep};
    d = N*Gamma_saved{iStep}*c + abs(N*G_total_saved{iStep})*ones(size(G_total_saved{iStep},2),1);
    R_assembled = polytope(N,d);
    % R_assembled should contain shrunk R
    R_shrunk = enlarge(R.timePoint.set{iStep},0.999);
    resPartial(end+1) = contains(R_assembled, R_shrunk);
    % R_assembled should not contain enlarged R
    R_enlarged = enlarge(R.timePoint.set{iStep},1.001);
    resPartial(end+1) = ~contains(R_assembled, R_enlarged);
end

%% compare reachable sets assembled from halfspace representation using G_tilde
% obtainlength of alpha vector
length_alpha = size(lin_A,2) - length(c);
for iStep = 1:length(R.timePoint.time)
    % assemble z
    z = [c;ones(length_alpha,1)];
    % assemble reachable set from computing linear constraints
    N = N_saved{iStep};
    d = [N*Gamma_saved{iStep}, G_tilde_saved{iStep}]*z;
    R_assembled = polytope(N,d);
    % R_assembled should contain shrunk R
    R_shrunk = enlarge(R.timePoint.set{iStep},0.999);
    resPartial(end+1) = contains(R_assembled, R_shrunk);
    % R_assembled should not contain enlarged R
    R_enlarged = enlarge(R.timePoint.set{iStep},1.001);
    resPartial(end+1) = ~contains(R_assembled, R_enlarged);
end

%% create alpha values: 
% small values slightly inside the unit box should fulfill the constaint,
% while values slightly outside the unit box should violate the constraint
% init alpha_small and alpha_large
alpha_small = zeros(length_alpha);
alpha_large = zeros(length_alpha);
% change only one value for alpha at a time
accuracy = 1e-12;
for i = 1:length_alpha
    % small values
    alpha_small(:,i) = ones(length_alpha,1) + accuracy;
    alpha_small(i,i) = 0.999;
    % large values
    alpha_large(:,i) = ones(length_alpha,1) + accuracy;
    alpha_large(i,i) = 1.001;
end

%% check constraints
for i = 1:length_alpha
    % small values should violate the constraint
    resPartial(end+1) = ~all(lin_A*[c;alpha_small(:,i)] <= lin_b);
    % large values should fulfill the constraint
    resPartial(end+1) = all(lin_A*[c;alpha_large(:,i)] <= lin_b);
end

%% final result
res = all(resPartial);

end


% Auxiliary functions -----------------------------------------------------

function [lin_A, lin_b, G_total_saved, Gamma_saved, N_saved, G_tilde_saved] = aux_linearConstraints(obj,union_y_a, params)

% generator matrices
G.X = generators(params.R0);
G.W = generators(params.W);
G.V = generators(params.V);
% numbers of generators
nrOfGen.X = size(G.X,2); % number of generators of X
nrOfGen.W = size(G.W,2); % number of generators of W
nrOfGen.V = size(G.V,2); % number of generators of V
% other variables
maxNrOfTimeSteps = ceil(params.tFinal/obj.dt); % maximum number of timeSteps
q = size(obj.C,1); % output dimension
n = size(obj.A,1); % system dimension

%% Compute constraints for finite time horizon
% Initialize values
Adi = eye(n);
lin_A = [];
lin_b = [];
Lambda = cell(maxNrOfTimeSteps+1,1); % further time step required to consider initial time
disturbanceMatrix = [];
Lambda{1} = obj.C;
LambdaSum = zeros(q,n);

% loop over all time steps
for k = 0:maxNrOfTimeSteps 
    % auxiliary values (for k>0); only for particular solution
    if k>0
        LambdaSum = LambdaSum + Lambda{k}; % C*[I + A_d + A_d^2 + ... + A_d^(k-1)]
        disturbanceMatrix(:,end+1:end+nrOfGen.W) = Lambda{k}*G.W; % C*[I, A_d, A_d^2, ..., A_d^(k-1)]*G.W
        % update i-th power of A
        Adi = obj.A*Adi; 
        % update Lambda
        Lambda{k+1} = obj.C*Adi; % Lambda_{k+1} is \Lambda_k in [1]
    else
        disturbanceMatrix = [];
    end
    % generator matrix G
    G_total = [Lambda{k+1}*G.X, disturbanceMatrix, G.V];
    % compute normal vectors of G
    N = aux_normalVectors(G_total);
    % nrOfConstraints
    nrOfConstraints = size(N,1);
    % Gamma
    Gamma = [Lambda{k+1}, LambdaSum, eye(q)];
    % compute G_tilde
    G_tilde = aux_compute_G_tilde(N, k, G, Lambda, nrOfGen);
    % partial A matrix
    A_partial = [-N*Gamma, -G_tilde];
    % partial b vector
    b_partial = -max(N*union_y_a{k+1}',[],2);
    % compose A and b for linear constraint by adding constraints of each
    % time step
    lin_A(end+1:end+nrOfConstraints,:) = A_partial; 
    lin_b(end+1:end+nrOfConstraints,1) = b_partial;
    
    % save intermediate values for additional checks
    G_total_saved{k+1} = G_total;
    Gamma_saved{k+1} = Gamma;
    N_saved{k+1} = N;
    G_tilde_saved{k+1} = G_tilde;
end
end

%% compute normal vectors of a zonotope with generator matrix G
% We do not use the halfspace conversion of zonotopes since the offset of
% the halfspaces to the origin is not required
function N = aux_normalVectors(G)
    % remove zero generators
    G = nonzeroFilter(G);
    % system dimension and number of generators
    [n,nrGen] = size(G);
    % get number of possible facets
    comb = combinator(nrGen,n-1,'c');
    % bypass bug in combinator (rows with all-zeros?!)
    comb = comb(any(comb,2),:);

    % build matrix of normal vectors
    N = zeros(length(comb(:,1)),n);
    for i=1:length(comb(:,1))
        indices=comb(i,:);
        Q=G(:,indices);
        v=ndimCross(Q);
        N(i,:)=v'/norm(v);
    end
    % second halg of normal vectors is symmetric
    N = [N; -N];
end

% computes auxiliary matrix
function G_tilde = aux_compute_G_tilde(M, k, G, Lambda, nrOfGen)
    % nr of constraints
    nrOfConstraints = size(M,1);
    %% A_X
    A_X = abs(M*Lambda{k+1}*G.X); % Lambda_{k+1} is \Lambda_k in [1]
    %% A_W
    if k > 0 % disturbance has only an effect for k>0
        % initialize result
        A_W = zeros(nrOfConstraints,nrOfGen.W);
        % loop over past time steps for input solution
        for i = 0:k-1
            % partial result
            A_W_partial = abs(M*Lambda{i+1}*G.W);
            % sum partial results
            A_W = A_W + A_W_partial;
        end
    else
        A_W = zeros(nrOfConstraints,nrOfGen.W);
    end
    %% A_V
    A_V = abs(M*G.V);
    %% assemble G_tilde
    G_tilde = [A_X, A_W, A_V];
end
        
% ------------------------------ END OF CODE ------------------------------
