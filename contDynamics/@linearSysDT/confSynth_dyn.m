function [params, union_y_a, cost] = confSynth_dyn(sys,params,options)
% confSynth_dyn - specific conformance synthesis method for linear 
%   discrete-time systems according to [1].
%
% Syntax:
%    params = confSynth_dyn(sys,params,options)
%
% Inputs:
%    sys - discrete-time linear system system
%    params - parameters defining the conformance problem
%    options - options for the conformance checking
%
% Outputs:
%    params - parameters solving the conformance problem
%    union_y_a - unified test cases (only deviation to nominal solution)
%    cost - cost of the optimization; required if system parameters are
%           optimized in an outer loop
%
% References:
%    [1] M. Althoff. Checking and Establishing Reachset Conformance in
%        CORA 2023. In Proc. of the 10th International Workshop on 
%        Applied Verification of Continuous and Hybrid Systems, 2023.
%    [2] Liu et al., "Guarantees for Real Robotic Systems: Unifying Formal
%        Controller Synthesis and Reachset-Conformant Identification", 202x.
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Stefan Liu
% Written:       29-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% Load variables

% generator matrices
G.X = generators(params.R0conf);
G.W = generators(params.W);
G.V = generators(params.V);
% numbers of generators
nrOfGen.X = size(G.X,2); % number of generators of X
nrOfGen.W = size(G.W,2); % number of generators of W
nrOfGen.V = size(G.V,2); % number of generators of V
% other variables
maxNrOfTimeSteps = ceil(params.tFinal/sys.dt); % maximum number of timeSteps
q = size(sys.C,1); % output dimension
n = size(sys.A,1); % system dimension
w = options.w; % weight vector

%% Unify all differences to the nominal solution; see y_a(k) in eq. 14 of [2]
[union_y_a,union_x_a] = conform_unifyTestCases(sys,params);

%% Corollary 1 in [2]: change time horizon to ensure conformance of an 
% infinite time horizon; requires that the entire state can be measured
if isinf(params.tFinal)
    % check enclosure of states after one time step in W
    
else
    %% Compute constraints for finite time horizon
    % Initialize values
    Adi = eye(n);
    lin_A = [];
    lin_b = [];
    Lambda = cell(maxNrOfTimeSteps+1,1); % further time step required to consider initial time
    disturbanceMatrix = [];
    Lambda{1} = sys.C;
    LambdaSum = zeros(q,n);

    % loop over all time steps
    for k = 0:maxNrOfTimeSteps 
        % auxiliary values (for k>0); only for particular solution
        if k>0
            LambdaSum = LambdaSum + Lambda{k}; % C*[I + A_d + A_d^2 + ... + A_d^(k-1)]
            disturbanceMatrix(:,end+1:end+nrOfGen.W) = Lambda{k}*G.W; % C*[I, A_d, A_d^2, ..., A_d^(k-1)]*G.W
            % update i-th power of A
            Adi = sys.A*Adi; 
            % update Lambda
            Lambda{k+1} = sys.C*Adi; % Lambda_{k+1} is \Lambda_k in [1]
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
    end
    % add constraints to enforce that alpha is positive
    length_alpha = nrOfGen.X + nrOfGen.W + nrOfGen.V;
    lin_A(end+1:end+length_alpha,:) = [zeros(length_alpha,2*n+q), -eye(length_alpha)];
    lin_b(end+1:end+length_alpha,1) = zeros(length_alpha, 1);
    % add some robustness
    rob = 1e-10;
    lin_b = lin_b - rob;
end

    
%% Optimization through linear programming or quadratic programming
switch options.norm
    %% Interval norm
    case 'interval'
        % Initialize norm
        iNorm = 0;
        % add norm values of further time steps
        for k = 0:maxNrOfTimeSteps
            % compute G_tilde
            G_tilde = aux_compute_G_tilde(eye(q), k, G, Lambda, nrOfGen);  
            % update norm
            iNorm = iNorm + w(k+1)*ones(1,q)*G_tilde;
        end
        % direction of optimization
        lin_f = [zeros(1,2*n+q), iNorm]';
        
        % Linear Programming
        opt = optimoptions('linprog','Display','none');
        [p_opt,cost,exitflag,output] = linprog(lin_f,lin_A,lin_b,[],[],[],[],opt);

    %% Frobenius norm
    case 'frob'
        % Initialize H_sum
        H_sum = 0;
        % add norm values of further time steps
        for k = 0:maxNrOfTimeSteps
            % compute H_tilde
            H_tilde = aux_compute_H_tilde(options.P, k, G, Lambda, nrOfGen);  
            % update norm
            H_sum = H_sum + w(k+1)*H_tilde;
        end
        % final H
        H = blkdiag(diag(zeros(1,2*n+q)), H_sum);
        
        % Quadratic Programming
        opt = optimoptions('quadprog','Display','none');
        [p_opt,cost,exitflag,output] = quadprog(H,[],lin_A,lin_b,[],[],[],[],[],opt);
        
    otherwise
        error('not supported')
end

% If no solution
if exitflag<=0
    error(output.message);
end

% Result
% p_opt is a vector of c_X, c_W, c_V, alpha_X, alpha_W, and alpha_V
% length of combined centers
dim_c = 2*n+q;
% assemble centers
c_X = p_opt(1:n);
c_W = p_opt(n+1:2*n);
c_V = p_opt(2*n+1:2*n+q);
% assemble alphas
alpha_X = p_opt(dim_c+1:dim_c+nrOfGen.X);
alpha_W = p_opt(dim_c+nrOfGen.X+1:dim_c+nrOfGen.X+nrOfGen.W);
alpha_V = p_opt(dim_c+nrOfGen.X+nrOfGen.W+1:end);

% new params for R0, W, and V
params.R0conf = zonotope([c_X, generators(params.R0conf)*diag(alpha_X)]);
params.W = zonotope([c_W, generators(params.W)*diag(alpha_W)]);
params.V = zonotope([c_V, generators(params.V)*diag(alpha_V)]);


end


% Auxiliary functions -----------------------------------------------------

%% compute normal vectors of a zonotope with generator matrix G
% We do not use the halfspace conversion of zonotopes since the offset of
% the halfspaces to the origin is not required
function N = aux_normalVectors(G)
    % remove zero generators
    G = nonzeroFilter(G);
    % system dimension and number of generators
    [n,nrGen] = size(G);
    % set is at least two-dimensional
    if n>1
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
        % second half of normal vectors is symmetric
        N = [N; -N];
    else
        % set "normal vectors" for one-dimensional interval
        N = [1; -1];
    end
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


% computes auxiliary matrix
function H_tilde = aux_compute_H_tilde(P, k, G, Lambda, nrOfGen)
    % dimension
    dim = size(G.X,1);
    %% H_X
    if ~isempty(G.X)
        H_X_aux = Lambda{k+1}*G.X;
        H_X = diag(diag(H_X_aux.'*P*H_X_aux));   % Lambda_{k+1} is \Lambda_k in [1]
    else
        H_X = [];
    end
    %% H_W
    if k > 0 % disturbance has only an effect for k>0
        % initialize inner result
        H_W_inner = zeros(dim,dim);
        % loop over past time steps for input solution
        for i = 0:k-1
            % partial result
            H_W_partial = Lambda{i+1}.'*P*Lambda{i+1};
            % sum partial results
            H_W_inner = H_W_inner + H_W_partial;
        end
        % H_W
        H_W = diag(diag(G.W.'*H_W_inner*G.W)); 
    else
        H_W = zeros(nrOfGen.W,nrOfGen.W);
    end
    %% A_V
    H_V = diag(diag(G.V.'*P*G.V));
    %% assemble H_tilde
    H_tilde = blkdiag(H_X, H_W, H_V);
end


% ------------------------------ END OF CODE ------------------------------
