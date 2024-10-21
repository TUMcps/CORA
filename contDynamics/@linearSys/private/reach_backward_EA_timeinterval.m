function timeInt = reach_backward_EA_timeinterval(linsys,params,options)
% reach_backward_EA_timeinterval - computes an inner approximation of the
%    time-interval maximal backward reachable set defined as
%       R(-tau) := { x_0 in R^n | forall u() exists w() exists t in tau:
%                                 x(t; x_0, u(), w()) in Xend }
%    via [1, Sec. VI-A-2)]
%
% Syntax:
%    timeInt = reach_backward_EA_timeinterval(linsys,params,options)
%
% Inputs:
%    linsys - linearSys object
%    params - model parameters
%    options - reachability settings
%
% Outputs:
%    timeInt - struct of time-interval solution
%
% References:
%    [1] M. Wetzlinger and M. Althoff, "Backward Reachability Analysis with
%        State Constraints for Linear Systems using Set Propagation", 2023.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachBackward

% Authors:       Mark Wetzlinger
% Written:       12-July-2023
% Last update:   11-August-2023 (MW, skip linear programs)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% time period and number of steps
tVec = params.tStart:options.timeStep:params.tFinal;
steps = length(tVec) - 1;
% time step size (shorter variable name)
dt = options.timeStep;

% initialize exponential matrix
expmat = exponentialMatrix(linsys.A);
% fast access/intuitive variable for state dimension
n = linsys.nrOfStates;

% initialize struct for backward reachable set
% array of time intervals
timeIntArray = cellfun(@(x) interval(x(1),x(2)), ...
    num2cell([-tVec(2:end)',-tVec(1:end-1)'],2),'UniformOutput',false);
timeInt.set = cell(steps,1);
timeInt.time = timeIntArray;

% is Rend an interval?
Rend_is_interval = isa(params.R0,'interval');

% number of inequalities in target set
if Rend_is_interval
    nrIneq = 2*n;
else
    nrIneq = size(params.R0.A,1);
end

% pre-compute particular solutions until start of considered time interval
% ...use 'closest' time step size to the one provided that partitions the
%    interval [0,params.tStart] into a whole number of steps
if ~withinTol(params.tStart,0)
    % integer number of steps from 0 to tStart
    steps_start = ceil(params.tStart/options.timeStep);
    dt_ = params.tStart/steps_start;

    % compute propagation matrix for one step
    expmat = init_eADeltat(expmat,dt_);
    % initialize cell arrays for propagation matrices
    expmat = init_eAtk(expmat,steps_start+steps);    

    % init directions for support function evaluation
    if Rend_is_interval
        l1_0 = [eye(n), -eye(n)];
        l2_0 = (l1_0' * expmat.eADeltat_back)';
    else
        l1_0 = params.R0.A';
        l2_0 = (params.R0.A * expmat.eADeltat_back)';
    end

    % compute inner- and outer-approximations of particular solution [0,dt_]
    [ZW_outer_0_dt,expmat] = oneStepPartSol(params.W,expmat,dt_,'outer');
    c_ZW = center(ZW_outer_0_dt); G_ZW = generators(ZW_outer_0_dt);
    [ZU_inner_0_dt,expmat] = oneStepPartSol(params.U,expmat,dt_,'inner');

    % propagate support function of particular solution due to disturbance
    % in 2*n axis-aligned directions
    ZW_outer_0_t_sF_l1 = zeros(nrIneq,1+steps_start+steps);
    ZW_outer_0_t_sF_l2 = zeros(nrIneq,1+steps_start+steps);
    for k=1:steps_start
        expmat = next_eAtk(expmat,k,'transpose');
        % propagation from start is numerically more robust
        l1 = expmat.eAtk_T{k} * l1_0;
        l2 = expmat.eAtk_T{k} * l2_0;
        % propagate support function value for first set of directions
        % (corresponding to P1)
        ZW_outer_0_t_sF_l1(:,k+1) = ZW_outer_0_t_sF_l1(:,k) ...
            + l1'*c_ZW + sum(abs(l1' * G_ZW),2);
        % propagate support function value for second set of directions
        % (corresponding to e^ADeltat * P1)
        ZW_outer_0_t_sF_l2(:,k+1) = ZW_outer_0_t_sF_l2(:,k) ...
            + l2'*c_ZW + sum(abs(l2' * G_ZW),2);
    end

    % propagate explicit particular solution
    ZU_inner_0_t = cell(1+steps_start+steps,1);
    ZU_inner_0_t{1} = zonotope(zeros(n,1));
    for k=1:steps_start
        expmat = next_eAtk(expmat,k);
        ZU_inner_0_t{k+1} = ZU_inner_0_t{k} + expmat.eAtk{k} * ZU_inner_0_dt;
    end

else
    steps_start = 0;
    % init directions and support function value of particular solution due
    % to distubances
    c_ZW_outer_0_t_sF_l1 = zeros(nrIneq,1+steps);
    c_ZW_outer_0_t_sF_l2 = zeros(nrIneq,1+steps);
    G_ZW_outer_0_t_sF_l1 = zeros(nrIneq,1+steps);
    G_ZW_outer_0_t_sF_l2 = zeros(nrIneq,1+steps);

    % init particular solution due to control inputs
    ZU_inner_0_t = cell(1+steps,1);
    ZU_inner_0_t{1} = zonotope(zeros(n,1));

    % init cell length for propagation matrices
    expmat = init_eAtk(expmat,steps);
end


if steps_start == 0 || ~withinTol(dt_,dt)
    % (re-)compute propagation matrix over one time step
    expmat = init_eADeltat(expmat,dt);

    % init directions for support function evaluation
    l1_0 = params.R0.A';
    l2_0 = (params.R0.A * expmat.eADeltat_back)';

    % re-compute particular solutions over [0,dt]
    [ZW_outer_0_dt,expmat] = oneStepPartSol(params.W,expmat,dt,'outer');
    c_ZW = center(ZW_outer_0_dt); G_ZW = generators(ZW_outer_0_dt);
    [ZU_inner_0_dt,expmat] = oneStepPartSol(params.U,expmat,dt,'inner');
end

% compute linear combination error
Rend_box_zono = zonotope(interval(params.R0));
% note: boxes always have n generators
mu = sqrt(n) * norm((expmat.eADeltat - eye(n))*generators(Rend_box_zono),2);

% compute curvature error
expmat = computeCurvatureErrors(expmat,dt);
C = expmat.F * Rend_box_zono;

% compute auxiliary polytopes
%    P1 = Rend - C - B(mu)
%    P2 = e^ADeltat*Rend - C - B(mu)
% where B(mu) is a ball of radius mu
P1 = minkDiff(polytope(params.R0),C);
P1 = polytope(P1.A,P1.b-mu);
% if representsa(P1,'emptySet')
%     keyboard
% end
P2 = minkDiff(expmat.eADeltat * polytope(params.R0),C);
P2 = polytope(P2.A,P2.b-mu);

% pre-compute parallelotopic outer-approximation of auxiliary polytopes
c1_box = center(Rend_box_zono);
G1_box = generators(Rend_box_zono);
c2_box = expmat.eADeltat * c1_box;
G2_box = expmat.eADeltat * G1_box;

% pre-compute support function in negative normal vectors of polytopes
P1_sF_minus = zeros(nrIneq,1);
P2_sF_minus = zeros(nrIneq,1);
if Rend_is_interval
    % we assume that the constraints (i=1,...,2n) are ordered such that
    %     A(i,:) = -A(n+i,:)
    P1_sF_minus = [P1.b(n+1:end); P1.b(1:n)];
    P2_sF_minus = [P2.b(n+1:end); P2.b(1:n)];
else
    for a=1:nrIneq
        P1_sF_minus(a) = supportFunc_(P1,-P1.A(a,:)','upper');
        P2_sF_minus(a) = supportFunc_(P2,-P2.A(a,:)','upper');
    end
end

% init directions
l1 = l1_0; l2 = l2_0;

% loop over all steps
for k=1:steps
    
    % propagate support function value of particular solution due to
    % disturbance in directions of P1 and P2
    c_ZW_outer_0_t_sF_l1(:,steps_start+k+1) = ...
        c_ZW_outer_0_t_sF_l1(:,steps_start+k) + l1'*c_ZW;
    G_ZW_outer_0_t_sF_l1(:,steps_start+k+1) = ...
        G_ZW_outer_0_t_sF_l1(:,steps_start+k) + sum(abs(l1' * G_ZW),2);
    c_ZW_outer_0_t_sF_l2(:,steps_start+k+1) = ...
        c_ZW_outer_0_t_sF_l2(:,steps_start+k) + l2'*c_ZW;
    G_ZW_outer_0_t_sF_l2(:,steps_start+k+1) = ...
        G_ZW_outer_0_t_sF_l2(:,steps_start+k) + sum(abs(l2' * G_ZW),2);

    % back-propagate directions for next iteration
    expmat = next_eAtk(expmat,steps_start+k+1,'transpose');
    l1 = expmat.eAtk_T{steps_start+k+1} * l1_0;
    l2 = expmat.eAtk_T{steps_start+k+1} * l2_0;

    % propagate propagation matrix
    expmat = next_eAtk(expmat,steps_start+k+1);

    % compute auxiliary constrained zonotope
    cZ_start = aux_minkDiffAndConZono(P1,c1_box,G1_box,P1_sF_minus,...
        c_ZW_outer_0_t_sF_l1(:,steps_start+k+1),...
        G_ZW_outer_0_t_sF_l1(:,steps_start+k+1),l1_0');
    cZ_end = aux_minkDiffAndConZono(P2,c2_box,G2_box,P2_sF_minus,...
        c_ZW_outer_0_t_sF_l2(:,steps_start+k+1),...
        G_ZW_outer_0_t_sF_l2(:,steps_start+k+1),l2_0');
    
    % propagate explicit particular solution due to control inputs
    if k < steps
        % skip last entry
        ZU_inner_0_t{steps_start+k+1,1} = ZU_inner_0_t{steps_start+k} ...
            + expmat.eAtk{steps_start+k} * ZU_inner_0_dt;
    end

    % time-interval solution (only save non-empty sets)
    if ~isemptyobject(cZ_start) && ~isemptyobject(cZ_end)
        % back-propagate propagation matrix
        expmat = next_eAtk(expmat,steps_start+k+1,'back');
        % compute time-interval solution
        CH_U = convHull(cZ_start,cZ_end) + (-ZU_inner_0_t{steps_start+k});
        timeInt.set{k,1} = expmat.eAtk_back{steps_start+k+1} * CH_U;
    else
        % once sets become empty, they are always empty
        % better: empty sets in reachSet object to show that sets are empty?
        break
    end

end

end


% Auxiliary functions -----------------------------------------------------

function cZ = aux_minkDiffAndConZono(...
    P,c_box,G_box,P_sF_minus,Z_sF_plus_c,Z_sF_plus_G,A_all)
% helper function to convert the Minkowski difference of a polytope P and a
% zonotope Z to a constrained zonotope cZ; this function exploits that we
% only need *some* zonotopic outer-approximation of the polytope; since P
% is constant over the loop in the main function, we pre-compute box(P) as
% a zonotope and use the corresponding zonotope center and generator matrix
% as input arguments; similarly, we use the support function values of Z
% in the directions of P as input arguments, which is considerably faster
% than computing the explicit particular solution
%
% References:
%    [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%        estimation and fault detection"
%    and polytope/conZonotope.m

    % number of constraints of the polytope
    nrIneq = size(A_all,1);

    % 1. compute Minkowski difference
    b_all = P.b - (Z_sF_plus_c + Z_sF_plus_G);
    % ... resulting polytope is <A_all,b_all> = minkDiff(P,Z)

    % 2. lower bound in the direction of i-th halfspace computed using
    %    sF(P-Z,l,'lower')
    %    = -sF(P-Z,-l,'upper')
    %    = -( sF(P,-l,'upper') - sF(Z,-l,'upper') )
    sigma = -(P_sF_minus - (Z_sF_plus_G - Z_sF_plus_c));
    
    % 3. construct constrained zonotope object according to eq. (21) in [1]
    G_ = [G_box, zeros(size(G_box,1),nrIneq)];
    A_ = [A_all*G_box, diag((b_all-sigma)./2)];
    b_ = (b_all+sigma)./2 - A_all*c_box;
    cZ = conZonotope(c_box,G_,A_,b_);
end

% ------------------------------ END OF CODE ------------------------------
