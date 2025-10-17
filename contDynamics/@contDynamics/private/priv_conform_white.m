function [params_new,fvals,p_opts,union_y_a,idzActive] = priv_conform_white(sys,params,options)
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
%           currently the interval norm ("interval", default) and the  
%           Frobenius norm ("frob") are supported
%       .cs.numRotations  - number of random rotations for the cost 
%           function "interval", default is zero
%       .cs.constraints - type of containment constraints, currently 
%           generator constraints ("gen", default) and halfspace 
%           constraints ("half") are supported
%       .cs.task - prediction task, currently regression ("reg", default) 
%           or classification ("class") are supported
%       .cs.w - weighting vector for the norm at different times for 
%           reachset conformance synthesis, the default value is a vector of ones
%       .cs.P - weighting matrix for the Frobenius norm, the default value 
%           is the identity matrix
%       .cs.robustness - robustness value for enforcing the 
%           containment constraints, the default value is $1e-9$.
%       .cs.updateDeriv - boolean, which can enforce the Jacobian 
%           and Hessian recomputation (required for linearization) at each 
%           nonlinear programming iteration
%       .cs.outMethod - method for outlier detection
%       .cs.numOutlier - number of data points allowed outside of
%           prediction set
%       .cs.idzOutlier - indizes of ouliers (data points that will be
%           considered in the cost computation but are allowed to violate
%           the constraints)
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
%    m_active - indizes of active data points at the final solution
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
        (~isprop(sys,'out_jacobian') || ~contains(func2str(sys.out_jacobian), "JacobianNN"))
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
n_y = sys.nrOfOutputs;

if isa(sys,'nonlinearSysDT') || isa(sys,'nonlinearARX') || contains(options.cs.recMethod,"rec")
    n_c = 0;
end

% set batch size and number of  recursions
if contains(options.cs.recMethod,"rec") % recursive identification
    assert(options.cs.cost == "interval") % TO-DO: remove? (should work?)
    if strcmp(options.cs.recMethod,"recUnSiCo")
        assert(options.cs.constraints == "half")
    end
    batchSize = options.cs.batchSize;
    numRecursions = ceil(length(testSuite)/batchSize);
else 
    batchSize = length(testSuite);
    numRecursions = 1;
end

% initalize variables for recursive identification
con.A_ineq = []; 
con.b_ineq = [];
lb_a = []; % initialize lower bound for alpha (will be set to default)
cost = 0;
fvals = zeros(1,numRecursions);

% update indizes and number of outliers
idzOutlier = options.cs.idzOutlier;
options.cs.numOutlier = options.cs.numOutlier - length(idzOutlier); % will be removed additionally

for batch=1:numRecursions
    idzBatch = (batch-1)*batchSize+1 : min(batch*batchSize,length(testSuite));

    % scale cost from previous batches with forgetting factor
    cost = options.cs.forgettingCost* cost;

    % compute the cost and constraint matrices for the current batch
    options.cs.idzOutlier = intersect(idzOutlier, idzBatch); % update outlier indizes
    [costBatch,conBatch,union_y_a] = aux_computeCostConstraints(sys,testSuite(idzBatch),G_X0,G_U,c_X0,c_U,options.cs,n_k,n_c);
    
    % combine cost and conastraint matrices with previous batches
    cost = cost + costBatch;
    con.A_ineq = [con.A_ineq; conBatch.A_ineq];
    con.b_ineq = [con.b_ineq; conBatch.b_ineq];
    con.A_eq = conBatch.A_eq;
    con.b_eq = conBatch.b_eq;

    % solve optimization problem
    n_m = length(idzBatch);
    if batch == 1
        % initialize optimal solution array
        p_opts = zeros(size(con.A_ineq,2),numRecursions);
    end
    [p_opts(:,batch),fvals(batch),problem] = aux_optimize(cost,con,options.cs,n_a,n_c,n_m,n_k,n_y,lb_a);

    if batch==ceil(length(testSuite)/batchSize)
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

if options.cs.determineActive && fvals(end) < inf
    % determine active data points
    if n_k > 1
        throw(CORAerror('CORA:specialError', ['Determining active measurements ' ...
            'not implemented for n_k > 1. Please use options.cs.determineActive = false.']))
    elseif isa(sys,'linearSysDT')
        throw(CORAerror('CORA:specialError', ['Determining active measurements ' ...
            'not implemented for linearSysDT. Please use nonlinearSysDT.']))
    end
    testSuite(idzOutlier) = [];
    idzActive = aux_determineActiveMeas(sys,p_opts(:,end),options,problem,testSuite,n_a,n_c); %,p_GO, G_X0, G_U, union_y_a,cost);
else
    idzActive = NaN;
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

function [p_opt,fval,problem,exitflag,output] = aux_optimize(cost,con,options_cs,n_a,n_c,n_m,n_k,n_y,lb_a)
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
    % quadratic program

    % outlier detection using frobenius norm cost not implemented
    assert(options_cs.numOutlier == 0)

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

elseif options_cs.numOutlier > 0 && strcmp(options_cs.outMethod,'MILP')
    % mixed-integer linear program
    % to allow outlier removal 
    if strcmp(options_cs.constraints,"half")
        throw(CORAerror('CORA:specialError', ['MILP outlier detection ' ...
            'not implemented for halfspace constraints. ' ...
            'Please use options.cs.constraints = ''gen''.']))
    end
    [p_opt,fval,exitflag,output] = aux_solveMILP(con,cost,lb,ub,options_cs,n_m,n_b+n_c,n_y);
    problem = NaN; % not returned
else  
    % linear program
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
end
if exitflag<=0
    % if no solution -> set fval and scaling factors to inf
    if options_cs.verbose
        output.message
    end
    disp("!!! No valid parameters found!!!")
    fval = inf;
    p_opt = [inf*ones(n_a,1); zeros(n_c+n_b,1)];
else
    % extract optimized scaling factors and center vectors
    if options_cs.verbose && output.constrviolation > 0
        CORAwarning('CORA:contDynamics',...
            "Small violation of the constraints during conformance identification.")
    end
    p_opt = p_opt(1:n_a+n_b+n_c);
end
end


function [p_opt,fval,exitflag,output] = aux_solveMILP(con,cost_alpha,lb,ub,options_cs,n_m,n_bc,n_y)
num_out = options_cs.numOutlier;
if strcmp(options_cs.task,"class")
    A11 = con.A_ineq(1:n_m*n_y,:);
    A12 = -diag(con.b_ineq(1:n_m*n_y))*kron(eye(n_m),ones(n_y,1));
    A21 = con.A_ineq(n_m*n_y+1:end,:);
    b2 = con.b_ineq(n_m*n_y+1:end);
    Aineq = [A11 A12; ...
        A21 zeros(size(A21,1),n_m); ...
        zeros(1,size(A11,2)) -ones(1,n_m)];
    bineq = [zeros(size(A11,1),1); ...
        b2; ...
        -(n_m-num_out)];
    Aeq = [];
    beq = [];
else
    %num_int = size(real(double(con.b_eq)),1);
    Aineq = [real(double(con.A_ineq)) zeros(size(con.A_ineq,1),n_m); ...
        zeros(1,size(con.A_ineq,2)) -ones(1,n_m)];
    bineq = [real(double(con.b_ineq)); -(n_m-num_out)];
    A12 = -diag(real(double(con.b_eq)))*kron(eye(n_m),ones(n_y,1));
    Aeq = [real(double(con.A_eq)) A12];
    beq = zeros(size(con.A_eq,1),1);
end
f = [double(cost_alpha) zeros(1, n_bc) zeros(1,n_m)];
intcon = length(f)-[n_m-1:-1:0];
lb = [lb; zeros(n_m,1)];
ub = [ub; ones(n_m,1)];

options = optimoptions('intlinprog','display','off'); %,'Algorithm','legacy');
[p_opt,fval,exitflag,output] = intlinprog(f',intcon,Aineq,bineq,Aeq,beq,lb,ub,[],options);
end


function idzActive = aux_determineActiveMeas(sys,p_opt,options,problem,testSuite,n_a,n_c) %,p_GO, G_X0, G_U, union_y_a,cost_alpha)
% solve optimization problem again

n_m = length(testSuite);
if all(p_opt == 0)
    fprintf(' (all predictions correct) \n')
    idzActive = [];
    return
end
n_b = size(problem.Aineq,2) - n_a - n_c;

if strcmp(options.cs.constraints,"gen")

    % run linear program to find active constraints
    problem_active = problem;
    alpha = p_opt(1:n_a,:);
    idz_akeep = alpha > 1e-6;
    idz_bkeep = repmat(idz_akeep,n_b/n_a,1);

    problem_active.f = [zeros(1,sum(idz_bkeep)) -ones(1, n_m)];
    problem_active.lb = [-repmat(alpha(idz_akeep),n_b/n_a,1); zeros(n_m,1)];
    problem_active.ub = [repmat(alpha(idz_akeep),n_b/n_a,1); repmat(min(alpha(idz_akeep)),n_m,1)];

    if strcmp(options.cs.task,"class")
        % classification tasks
        Aeq = problem.Aineq(1:sys.nrOfOutputs*n_m,:);
        beq = problem.bineq(1:sys.nrOfOutputs*n_m);
        Aineq = problem.Aineq(sys.nrOfOutputs*n_m+1:end,:);
        bineq = problem.bineq(sys.nrOfOutputs*n_m+1:end);
    else
        % regression tasks
        Aeq = problem.Aeq;
        beq = problem.beq;
        Aineq = problem.Aineq;
        bineq = problem.bineq;
    end
    beq = beq - Aeq(:,1:n_a)*alpha  - Aeq(:,n_a+1:n_a+n_c)*p_opt(n_a+1:n_a+n_c);
    Aeq = Aeq(:,n_a+n_c+1:end);
    Aeq = [Aeq(:,idz_bkeep) zeros(size(Aeq,1),n_m)];

    bineq = -Aineq(:,1:n_a)*alpha;
    idx_bzero = abs(bineq) <= 1e-6;
    bineq(idx_bzero) = [];
    Aineq = Aineq(:,n_a+1:end);
    f = double(idz_akeep);
    Aineq = [Aineq(:,idz_bkeep) repmat(kron(eye(n_m),f),2,1)];
    Aineq(idx_bzero,:) = [];
    if strcmp(options.cs.task,"class")
        % classification tasks
        problem_active.Aeq = [];
        problem_active.beq = [];
        problem_active.Aineq = [Aeq; Aineq];
        problem_active.bineq = [beq; bineq];
    else
        % regression tasks
        problem_active.Aeq = Aeq;
        problem_active.beq = beq;
        problem_active.Aineq = Aineq;
        problem_active.bineq = bineq;
    end

    % solve
    [p_optCon,fvalCon,exitflag,output] = CORAlinprog(problem_active);
    try
        boolActive = p_optCon(end-n_m+1:end) < 1e-6;
    catch
        fprintf(" (optimization error: %s)", output.message)
        boolActive = ones(n_m,1);
    end
    if sum(boolActive)==0
        fprintf(" (no active constraints found)")
        boolActive = ones(n_m,1);
    end

else
    % find all active testcases
    activeC = (problem.Aineq*p_opt(1:n_a+n_b+n_c,:)-problem.bineq >= -options.cs.robustness);
    if contains(options.cs.cost,"class")
        num_constr2 = n_m * sys.nrOfOutputs;
        activeC_m1 = reshape(activeC(1:end-num_constr2),[],n_m);
        activeC_m2 = reshape(activeC(end-num_constr2+1:end),[],n_m);
        boolActive = any(activeC_m1,1) & any(activeC_m2,1);
    else        
        activeC_m = reshape(activeC, [], n_m);
        boolActive = (sum(activeC_m,1)> 0)';
    end
end
idzActive = find(boolActive > 0);

% Check new solution only using the active constraints
if false
    con = priv_conform_whiteCon(sys, testSuite(boolActive), p_GO(boolActive), G_X0, G_U, options.cs, union_y_a(boolActive));
    % solve optimization problem
    % dimensions of the optimization variables alpha, beta and center
    n_b = size(con.A_ineq,2) - n_a - n_c;
    % optimization bounds
    lim_b = inf;
    lb = [options.cs.a_min*ones(n_a,1); ...
        -options.cs.cp_lim*ones(n_c,1); ...
        -lim_b*ones(n_b,1)];
    ub = [options.cs.a_max*ones(n_a,1); ...
        options.cs.cp_lim*ones(n_c,1); ...
        lim_b*ones(n_b,1)];
    cost = [double(cost_alpha) zeros(1, n_c) zeros(1, n_b)];

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
    [p_opt2,fval2,exitflag,output] = CORAlinprog(problem);
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

if ~isempty(options_cs.idzOutlier)
    % remove outliers specified by idzOutlier
    testSuite(options_cs.idzOutlier) = [];
    if iscell(p_GO)
        p_GO(options_cs.idzOutlier) = [];
    else
        p_GO.y(:,:,options_cs.idzOutlier) = [];
    end
    if iscell(union_y_a)
        union_y_a(options_cs.idzOutlier) = [];
    else
        union_y_a(:,:,options_cs.idzOutlier) = [];
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

elseif strcmp(options_cs.cost,"interval") % (rotated) interval norm
    cost_rot = 0;
    R = [];
    n = size(p_GO.C{1},1);
    tol=1e-6;

    % compute rotation matrices and stack vertically
    for j=1:options_cs.numRotations
        M = zeros(n); % prealloc
        % gram-schmidt on random column vectors
        vi = randn(RandStream('mt19937ar'),n,1);
        % The n-dimensional normal distribution has spherical symmetry.
        % Thus, after normalization the drawn vectors would be uniformly
        % distributed on the n-dimensional unit sphere.
        M(:,1) = vi ./ norm(vi);

        for i=2:n
            nrm = 0;
            while nrm<tol
                vi = randn(n,1);
                vi = vi -  M(:,1:i-1)  * ( M(:,1:i-1).' * vi )  ;
                nrm = norm(vi);
            end
            M(:,i) = vi ./ nrm;
        end %i
        R = [R;M];
    end

    for k = 1:n_k
        % loop through time steps of the trajectory
        sum_U = 0;
        sum_rot = 0;
        for j = 1:k
            sum_U = sum_U + abs(p_GO.D{k,j} * G_U);
            if options_cs.numRotations > 0
                sum_rot = sum_rot + sum(abs(R*p_GO.D{k,j} * G_U),1);
            end
        end
        if options_cs.numRotations > 0
            cost_rot = cost_rot + [sum_rot sum(abs(R*p_GO.C{k} * G_X0),1)];
        end
        cost = cost + options_cs.w(k) * [abs(p_GO.C{k} * G_X0) sum_U];
    end
    cost = sum(cost,1) + cost_rot;
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
