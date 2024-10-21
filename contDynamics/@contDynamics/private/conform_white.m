function [params_new,fval,p_opt,union_y_a] = conform_white(sys,params,options)
% conform_white - specific conformance synthesis method for nonlinear 
%    discrete-time systems according to [1].
%
% Syntax:
%    [params_new,fval,p_opt,union_y_a] = conform_white(sys,params,options)
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
%       .cs.robustnessMargin - robustness value for enforcing the 
%           containment constraints, the default value is $1e-9$.
%       .cs.derivRecomputation - boolean, which can enforce the Jacobian 
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
%    union_y_a - n_m x 1 cell array with n_k x n_y x n_s output arrays
%       with union_y_a{m} = testSuite{m}.y - y_nom
%
% Example: 
%    -
% 
% References:
%    [1] L. Luetzow and M. Althoff, "Reachset-conformant System 
%        Identification," arXiv, 2024.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Laura Luetzow, Matthias Althoff
% Written:       25-July-2023
% Last update:   28-March-2024 (LL, restructuring)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if options.cs.derivRecomputation && ...
        (isa(sys, 'nonlinearSysDT') || isa(sys, 'nonlinearARX')) && ...
        ~contains(func2str(sys.mFile), "predict")
    % derivative recomputation since model is changed
    derivatives(sys,options);
end

n_k = ceil(round(params.tFinal/sys.dt,2)) + 1; % maximum number of timeSteps
testSuite = params.testSuite;
G_X0 = generators(params.R0);
G_U = generators(params.U);

% compute the linearized system matrices and the cost for each test case m
p_GO = cell(length(testSuite),1);
cost = 0;

union_y_a = cell(length(testSuite),1);
for m = 1 : length(testSuite)
    u_nom = testSuite{m}.u + center(params.U)';
    x0_nom = testSuite{m}.initialState + center(params.R0);
    p_GO{m} = computeGO(sys, x0_nom, permute(u_nom,[2,1,3]), n_k);
    union_y_a{m} = testSuite{m}.y - permute(p_GO{m}.y,[2 1 3]);

    % compute cost of Y_X0U
    cost = cost + aux_conformCost(p_GO{m}, G_X0, G_U, options.cs);
end

% compute the constraint parameters
con = conform_whiteCon(sys, testSuite, p_GO, G_X0, G_U, options.cs, union_y_a);

% solve optimization problem
% dimensions of the optimization variables alpha, beta and center
n_a = options.cs.n_a;
n_c = options.cs.n_c;
n_b = size(con.A_ineq,2) - n_a - n_c;

% optimization bounds
lim_b = inf;
lb = [options.cs.a_min*ones(n_a,1); ...
    -options.cs.cp_lim*ones(n_c,1); ...
    -lim_b*ones(n_b,1)];
ub = [options.cs.a_max*ones(n_a,1); ...
    options.cs.cp_lim*ones(n_c,1); ...
    lim_b*ones(n_b,1)];


% compute parameters with LP or QP using the linear constraints
con.b_ineq = con.b_ineq - options.cs.robustnessMargin;
if strcmp(options.cs.cost,"frob")
    % quadratic cost and linear constraints
    cost = blkdiag(double(cost), zeros(n_c+n_b, n_c+n_b));
    cost_lin = zeros(1, size(cost,1));

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
    p0 = [ones(n_a, 1); zeros(n_c+n_b, 1)];
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

elseif strcmp(options.cs.cost,"interval")  % interval norm
    % linear cost and constraints
    cost = [double(cost) zeros(1, n_c) zeros(1, n_b)];

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
    % If no solution
    if exitflag<=0
        if options.cs.verbose
            output.message
            disp("!!! No valid parameters found!!!")
        end
        fval = inf;
        p_opt = [zeros(n_c, 1); ones(n_a, 1)];
    else
        if options.cs.verbose && output.constrviolation > 0
            CORAwarning('CORA:contDynamics',...
                "Small violation of the constraints during conformance identification.")
        end
        p_opt = p_opt(1:n_a+n_c);
    end
end

% create the uncertainty sets from the estimated parameters p_opt
p_a = p_opt(1:n_a);
c_X0 = center(params.R0);
c_U = center(params.U);
if isa(sys,'linearSysDT') || isa(sys,'linearARX')
    p_c = p_opt(n_a+1:n_a+n_c);
    c_X0 = c_X0 + p_c(1:size(G_X0,1)); %?????
    c_U = c_U + p_c(size(G_X0,1)+1:end);
end

% build the new parameter object
G_X0 = generators(params.R0) * diag(p_a(1:size(G_X0,2)));
G_U = generators(params.U) * diag(p_a(size(G_X0,2)+1:end));
params_new = params;
params_new.U = zonotope(c_U,G_U);
params_new.R0 = zonotope(c_X0,G_X0);

end


% Auxiliary functions -----------------------------------------------------

function cost = aux_conformCost(p_GO,G_X0,G_U,options_cs)
% cost function for reachset conformance testing 

n_k = size(p_GO.y,2);
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
        cost = cost + options_cs.w(k) * blkdiag(diag_X0, diag_U);
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

% ------------------------------ END OF CODE ------------------------------
