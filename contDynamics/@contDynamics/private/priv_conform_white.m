function [params_new,fvals,p_opts,union_y_a] = priv_conform_white(sys,params,options)
% priv_conform_white - specific conformance synthesis method for 
%    discrete-time systems according to [1] and [2].
%
% Syntax:
%    [params_new,fval,p_opt,union_y_a] = priv_conform_white(sys,params,options)
%
% Inputs:
%    sys - discrete-time nonlinear system 
%    params - parameters defining the conformance problem
%    options - options for the conformance checking
%       .cs.cost - cost function for reachset conformance synthesis, 
%           currently the interval norm ("interval") and the Frobenius norm 
%           ("frob") are supported
%       .cs.w - weighting vector for the norm at different times for 
%           reachset conformance synthesis, the default value is a vector of ones
%       .cs.P - weighting matrix for the Frobenius norm, the default value 
%           is the identity matrix
%       .cs.constraints - type of containment constraints, currently 
%           halfspace constraints ("half") and generator constraints 
%           ("gen") are supported
%       .cs.robustness - robustness value for enforcing the 
%           containment constraints, the default value is $1e-9$.
%       .cs.updateDeriv - boolean, which can enforce the Jacobian 
%           and Hessian recomputation (required for linearization) at each 
%           nonlinear programming iteration
%       .cs.a_min - lower limit for the identified scaling factors, default 
%           is 0
%       .cs.a_max - upper limit for the identified scaling factors, default 
%           is infinite
%       .cs.cp_lim - upper limit for the absolute value of the identified 
%           center vectors and parameters, default is infinite
%       .cs.verbose - boolean, which can suppress displaying output
%
% Outputs:
%    params_new - parameters solving the conformance problem
%    fval - conformance cost
%    p_opt - estimated parameters containing the scaling factors alpha 
%       and the change of the center vectors Delta c with 
%       p_opt = [alpha_x' alpha_u' dc_x' dc_u']'
%    union_y_a - n_m x 1 cell array with dim_y x n_k x n_s output arrays
%       with union_y_a{m} = testSuite(m).y - y_nom
%
% Example: 
%    -
% 
% References:
%    [1] L. Luetzow and M. Althoff, "Reachset-conformant System 
%        Identification," arXiv, 2025.
%    [2] L. Luetzow and M. Althoff, "Recursive identification of reachset-
%        conformant models using constraint underapproximation," in Proc. 
%        of the American Control Conference, 2025.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Laura Luetzow, Matthias Althoff
% Written:       25-July-2023
% Last update:   28-March-2024 (LL, restructuring)
%                03-March-2025 (LL, add recursive identification)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if options.cs.updateDeriv && ...
        (isa(sys,'nonlinearSysDT') || isa(sys,'nonlinearARX')) && ...
        ~contains(func2str(sys.mFile),"predict")
    % derivative recomputation since model is changed
    derivatives(sys,options);
end

n_k = ceil(round(params.tFinal/sys.dt,2)) + 1; % maximum number of timeSteps
testSuite = params.testSuite;
G_X0 = generators(params.R0);
G_U = generators(params.U);
c_X0 = center(params.R0);
c_U = center(params.U);

% dimensions of the optimization variables
n_a = size(params.U.generators,2) + size(params.R0.generators,2);
n_c = size(params.U.generators,1) + size(params.R0.generators,1);

if isa(sys,'nonlinearSysDT') || isa(sys,'nonlinearARX') || contains(options.cs.recMethod,"rec")
    n_c = 0;
end

% set batch size and number of  recursions
if contains(options.cs.recMethod,"rec") % recursive identification
    assert(options.cs.cost == "interval") % TO-DO: remove? (should work?)
    if strcmp(options.cs.recMethod,"recUnSiCo")
        assert(options.cs.constraints == "half")
    end
    batch_size = options.cs.batchSize;
    num_recursions = ceil(length(testSuite)/batch_size);
else 
    batch_size = length(testSuite);
    num_recursions = 1;
end

% initalize variables for recursive identification
con.A_ineq = []; 
con.b_ineq = [];
lb_a = []; % initialize lower bound for alpha (will be set to default)
cost = 0;
fvals = zeros(1,num_recursions);
p_opts = zeros(n_a+n_c,num_recursions);

for batch=1:num_recursions
    batch_idz = (batch-1)*batch_size+1 : min(batch*batch_size,length(testSuite));

    % scale cost from previous batches with forgetting factor
    cost = options.cs.forgettingCost* cost;

    % compute the cost and constraint matrices for the current batch
    [cost_batch,con_batch,union_y_a] = aux_computeCostConstraints(sys,testSuite(batch_idz),G_X0,G_U,c_X0,c_U,options.cs,n_k,n_c);
    
    % combine cost and conastraint matrices with previous batches
    cost = cost + cost_batch;
    con.A_ineq = [con.A_ineq; con_batch.A_ineq];
    con.b_ineq = [con.b_ineq; con_batch.b_ineq];
    con.A_eq = con_batch.A_eq;
    con.b_eq = con_batch.b_eq;

    % solve optimization problem
    [p_opts(:,batch),fvals(batch),problem] = aux_optimize(cost,con,options.cs,n_a,n_c,lb_a);

    if batch==ceil(length(testSuite)/batch_size)
        % maximum number of iterations reached
        break
    end

    % update the constraints for the recursive methods
    if strcmp(options.cs.recMethod,"recUnSiCo")
        [con.A_ineq,con.b_ineq] = aux_UnSiCoId(p_opts(:,batch),con.A_ineq,con.b_ineq,cost,problem,options);

    else % "recAlpha"
        con.A_ineq = [];
        con.b_ineq = [];

        %set lower bound to estimated alpha
        lb_a = options.cs.forgetting*p_opts(1:n_a,batch);
    end
end

if n_c > 0
    % update the center vectors
    p_c = p_opts(n_a+1:n_a+n_c, end);
    c_X0 = c_X0 + p_c(1:size(G_X0,1));
    c_U = c_U + p_c(size(G_X0,1)+1:end);
end

% build the uncertainty sets from the optimized parameters
p_a = p_opts(1:n_a,end);
G_X0 = generators(params.R0) * diag(p_a(1:size(G_X0,2)));
G_U = generators(params.U) * diag(p_a(size(G_X0,2)+1:end));
params_new = params;
G_U(isnan(G_U)) = inf; % invalid results
G_X0(isnan(G_X0)) = inf; % invalid results
params_new.U = zonotope(c_U,G_U);
params_new.R0 = zonotope(c_X0,G_X0);

end


% Auxiliary functions -----------------------------------------------------

function [p_opt,fval,problem,exitflag,output] = aux_optimize(cost,con,options_cs,n_a,n_c,lb_a)
% create the otimization problem and solve it for the scaling factors and
% center vectors

% dimension of auxiliary optimization variable beta
n_b = size(con.A_ineq,2) - n_a - n_c;

% set lower and upper bound for the optimization variables
lim_b = inf;
if isempty(lb_a)
    lb_a = options_cs.a_min*ones(n_a,1);
end
lb = [lb_a; ...
    -options_cs.cp_lim*ones(n_c,1); ...
    -lim_b*ones(n_b,1)];
ub = [options_cs.a_max*ones(n_a,1); ...
    options_cs.cp_lim*ones(n_c,1); ...
    lim_b*ones(n_b,1)];

% compute parameters with LP or QP using the linear constraints
con.b_ineq = con.b_ineq - options_cs.robustness;
if strcmp(options_cs.cost,"frob")
    % quadratic cost and linear constraints
    cost = blkdiag(double(cost),zeros(n_c+n_b,n_c+n_b));
    cost_lin = zeros(1,size(cost,1));

    % construct problem
    problem = struct;
    problem.H = cost;
    problem.f = cost_lin;
    problem.Aineq = con.A_ineq;
    problem.bineq = con.b_ineq;
    problem.Aeq = con.A_eq;
    problem.beq = con.b_eq;
    problem.lb = lb;
    problem.ub = ub;

    % initial point
    p0 = [ones(n_a,1); zeros(n_c+n_b,1)];
    problem.x0 = p0;

    [p_opt,fval,exitflag,output] = CORAquadprog(problem);
    % If no solution
    if exitflag <= 0
        output.message;
        fval = inf;
        p_opt = p0;
    else
        p_opt = p_opt(1:n_a+n_c);
    end

elseif strcmp(options_cs.cost,"interval")  % interval norm
    % linear cost and constraints
    cost = [double(cost) zeros(1,n_c) zeros(1,n_b)];

    % make constraint matrices real and double
    problem = struct;
    problem.f = cost;
    problem.Aineq = real(double(con.A_ineq));
    problem.bineq = real(double(con.b_ineq));
    problem.Aeq = real(double(con.A_eq));
    problem.beq = real(double(con.b_eq));
    problem.lb = lb;
    problem.ub = ub;

    % call linprog
    [p_opt,fval,exitflag,output] = CORAlinprog(problem);
    if exitflag<=0
        % if no solution -> set fval and scaling factors to inf
        if options_cs.verbose
            output.message
        end
        disp("!!! No valid parameters found!!!")
        fval = inf;
        p_opt = [zeros(n_c,1); inf*ones(n_a,1)];
    else
        % extract optimizes scaling factors and center vectors
        if options_cs.verbose && output.constrviolation > 0
            CORAwarning('CORA:contDynamics',...
                "Small violation of the constraints during conformance identification.")
        end
        p_opt = p_opt(1:n_a+n_c);
    end
end
end

function [cost,con,union_y_a] = aux_computeCostConstraints(sys,testSuite,G_X0,G_U,c_X0,c_U,options_cs,n_k,n_c)
% compute the cost coefficient and the constraint matrices

p_GO = cell(length(testSuite),1);
cost = 0;
union_y_a = cell(length(testSuite),1);

% compute the linearized system matrices and the cost for each test case m
if isa(sys, 'linearARX') || isa(sys, 'linearSysDT')
    % more efficient computation for linear systems
    u_nom = NaN(size(testSuite(1).u,1),size(testSuite(1).u,2), length(testSuite));
    x0_nom = zeros(size(testSuite(1).x,1),1, length(testSuite));

    for m = 1 : length(testSuite)
        % compute GO model for each test case
        if testSuite(m).n_k > size(u_nom,2)
            % augment u_nom with NaN values if trajectory m is longer
            u_nom = [u_nom NaN(size(u_nom,1), testSuite(m).n_k-size(u_nom,2),size(u_nom,3))];
        end
        u_nom(:,1:testSuite(m).n_k,m) = testSuite(m).u + c_U;
        x0_nom(:,1,m) = testSuite(m).x(:,1,1) + c_X0;
    end
    p_GO = computeGO(sys,x0_nom,u_nom,n_k);
    union_y_a = zeros(size(p_GO.y,1),0,0);
    for m = 1 : length(testSuite)
        union_y_a_m = testSuite(m).y - p_GO.y(:,1:size(testSuite(m).y,2),m);
        % fill with NaN if unequal trajectory lengths
        if size(union_y_a,2) > size(union_y_a_m,2)
            union_y_a_m = [union_y_a_m NaN(size(union_y_a_m,1),size(union_y_a,2)-size(union_y_a_m,2),size(union_y_a_m,3))];
        elseif size(union_y_a,2) < size(union_y_a_m,2)
            union_y_a = [union_y_a NaN(size(union_y_a,1),size(union_y_a_m,2)-size(union_y_a,2),size(union_y_a,3))];
        end

        union_y_a = cat(3, union_y_a, union_y_a_m);
    end

    % compute cost of Y_X0U
    cost = cost + aux_conformCost(p_GO,G_X0,G_U,options_cs);
else
    for m = 1 : length(testSuite)
        % compute GO model for each test case
        u_nom = testSuite(m).u + c_U;
        x0_nom = testSuite(m).x(:,1,1) + c_X0;
        p_GO{m} = computeGO(sys,x0_nom,u_nom,n_k);
        union_y_a{m} = testSuite(m).y - p_GO{m}.y(:,1:size(testSuite(m).y,2));

        % compute cost of Y_X0U
        cost = cost + aux_conformCost(p_GO{m},G_X0,G_U,options_cs);
    end
end

% compute the constraint parameters
con = priv_conform_whiteCon(sys,testSuite,p_GO,G_X0,G_U,options_cs,union_y_a,n_c);
end

function cost = aux_conformCost(p_GO,G_X0,G_U,options_cs)
% cost function for reachset conformance testing

n_k = size(p_GO.C,1);
cost = 0;

if strcmp(options_cs.cost,"frob")  % frobenius norm
    for k = 1:n_k
        % loop through time steps of the trajectory
        sum_U = 0;
        for j = 1:k
            sum_U = sum_U + p_GO.D{k,j}' * options_cs.P * p_GO.D{k,j};
        end
        diag_U = diag(diag(G_U' * sum_U * G_U));
        diag_X0 = diag(diag(G_X0' *p_GO.C{k}' * options_cs.P * p_GO.C{k} * G_X0));
        cost = cost + options_cs.w(k) * blkdiag(diag_X0,diag_U);
    end

elseif strcmp(options_cs.cost,"interval")  % interval norm
    for k = 1:n_k
        % loop through time steps of the trajectory
        sum_U = 0;
        for j = 1:k
            sum_U = sum_U + abs(p_GO.D{k,j} * G_U); %--> G_U = [Gu Gw Gv]
        end
        cost = cost + options_cs.w(k) * [abs(p_GO.C{k} * G_X0) sum_U];
    end
    cost = sum(cost,1);
end

end

function [A,b] = aux_UnSiCoId(p_opt,A_new,b_new,cost,problem,options)
P_opt = repmat(p_opt,1,options.cs.numPoints);
for j=2:options.cs.numPoints
    problem.f = cost.*(2*rand(size(cost)));
    [p_opt2,~,exitflag,~] = CORAlinprog(problem);

    % If no solution
    if exitflag<=0
        p_opt2 = p_opt;
    end
    P_opt(:,j) = p_opt2;
end

% reduce number of cosntraints with clustering
% normalize
A_norm = A_new ./ -sum(A_new,2);
b_norm = b_new ./ -sum(A_new,2);

% recursive clutering
[A,b] = aux_UnSiCo(A_norm,b_norm,P_opt,options.cs.numCluster);
b = options.cs.forgetting*b;
end

function [A,b] = aux_UnSiCo(A_all,b_all,p_opt,num_cluster)
% underapproximation of similar constraints

if num_cluster >= size(A_all,1) / 2
    A = A_all;
    b = b_all;
    return
end
idx = kmeans([A_all b_all],num_cluster);

% overapproximate each cluster
A = [];
b = [];
for idx_i= 1:num_cluster
    A_cl_i = A_all(idx==idx_i,:);
    b_cl_i = b_all(idx==idx_i,:);

    if any(max(A_cl_i,[],1)*p_opt-min(b_cl_i,[],1) > 0,'all')
        [A_res,b_res] = aux_UnSiCo(A_cl_i,b_cl_i,p_opt,num_cluster);
        A = [A; A_res];
        b = [b; b_res];
    else
        A = [A; max(A_cl_i,[],1)];
        b = [b; min(b_cl_i,[],1)];
    end
end
end

% ------------------------------ END OF CODE ------------------------------
