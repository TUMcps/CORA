function [timeInt,timePoint,res] = priv_reach_krylov(linsys,params,options)
% priv_reach_krylov - computes the reachable set for linear systems using
% 	 the krylov reachability algorithm for linear systems [1]
%
% Syntax:
%    [timeInt,timePoint,res] = priv_reach_krylov(linsys,params,options)
%
% Inputs:
%    linsys - linearSys object
%    params - model parameters
%    options - options for the computation of reachable sets
%
% Outputs:
%    timeInt - array of time-interval output sets
%    timePoint - array of time-point output sets
%    res - true/false (only if specification given)
%
% Example:
%    -
%
% References:
%    [1] M. Althoff. "Reachability analysis of large linear systems with
%        uncertain inputs in the Krylov subspace", IEEE Transactions on
%        Automatic Control 65 (2), pp. 477-492, 2020.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Maximilian Perschl, Matthias Althoff, Mark Wetzlinger
% Written:       26-June-2019
% Last update:   19-November-2022 (MW, modularize specification check)
%                23-April-2025 (MP, refactor including new error bound)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% if a trajectory should be tracked
if isfield(params,'uTransVec')
    params.uTrans = params.uTransVec(:,1);
end

% log information
verboseLog(options.verbose,1,params.tStart,params.tStart,params.tFinal);

%initialize reachable set computations
[linsys,params,options] = priv_initReach_Krylov(linsys,params,options);

%time period
tVec = params.tStart:options.timeStep:params.tFinal;
steps = length(tVec);

%create options.t
options.t = params.tStart;

% initialize arguments for the output equation
timeInt.set = cell(steps-1,1);
timeInt.time = cell(steps-1,1);
timePoint.set = cell(steps,1);
timePoint.time = num2cell(tVec');

% fill in values for initial set
timePoint.set{1} = linsys.C*params.R0;
timePoint.time{1} = 0;


% loop over all reachability steps
for k=2:steps
    
    % increment time
    options.t = tVec(k);

    % propagate different sets using exponential matrix in Krylov subspace
    Rnext = aux_propagation(linsys,params,options);

    % calculate the output set
    timeInt.set{k-1} = Rnext.ti;
    timeInt.time{k-1} = interval(tVec(k-1),tVec(k));
    timePoint.set{k} = Rnext.tp;
    timePoint.time{k} = options.t;
    
    % safety property check
    if isfield(options,'specification')
        [res,timeInt,timePoint] = priv_checkSpecification(...
            options.specification,[],timeInt,timePoint,k-1);
        if ~res; return; end
    end
    
    % log information
    verboseLog(options.verbose,k,tVec(k),params.tStart,params.tFinal);
    
    % if a trajectory should be tracked
    if isfield(params,'uTransVec')
        params.uTrans = params.uTransVec(:,k);
        % input trajectory not yet implemented
    end
   
end

% specification fulfilled
res = true;
end


% Auxiliary functions -----------------------------------------------------

function Rnext = aux_propagation(sys,params,options)
% Propagate all set solutions 
%           Rnext - Reachable set of next timestep
%               R.tp - timepoint solution
%               R.ti - timeinterval solution

% retrieve reachable set
U_0_sol = sys.krylov.Rpar_proj;
Rhom_tp_prev = sys.krylov.Rhom_tp_prev;

% propagate
% homogeneous solution
[Htp,TIE_new_proj] = aux_propagate_HomSol(sys.krylov.state.c_sys_proj,...
                                sys.krylov.state.g_sys_proj,options.timeStep,options.t - options.timeStep,params.R0,options.taylorTerms,options.krylovError,sys.nrOfOutputs);

% uncertain input solution
if ~isnumeric(U_0_sol)
    [RV_proj,~] = priv_exponential_Krylov_projected_linSysInput(sys.krylov.input.c_sys_proj,sys.krylov.input.g_sys_proj,sys.krylov.RV_0,sys.nrOfOutputs,options);
else
    % for the first step, we already computed the input solution in
    % init_reach_krylov
    RV_proj = sys.krylov.Rpar_proj_0;
end
% constant input solution
% dummy R0 to get norm correct
[RTrans_new_proj,~] = ...
    priv_exponential_Krylov_projected_linSysInput(sys.krylov.uTrans_sys,[],zonotope(1),sys.nrOfOutputs,options);


%% compute new solutions for this timestep

% we save U_0_sol as interval because otherwise we would have way too many
% generators to justify the loss in precision
U_0_sol = U_0_sol + interval(RV_proj);
U_0_sol = reduce(U_0_sol,options.reductionTechnique,options.zonotopeOrder);


% next time step solution
Rhom_tp_proj = Htp + RTrans_new_proj;
R_tp_proj = Rhom_tp_proj + zonotope(U_0_sol);
% affine time interval solution
TI_apx = enclose(Rhom_tp_prev,Rhom_tp_proj);
% complete time-interval solution
R_ti_proj = TI_apx + U_0_sol + TIE_new_proj + sys.krylov.inputCorr;

% order reduction
Rnext.ti = reduce(R_ti_proj,options.reductionTechnique,options.zonotopeOrder);
Rnext.tp = reduce(R_tp_proj,options.reductionTechnique,options.zonotopeOrder);

% save result
sys.krylov.Rhom_tp_prev = Rhom_tp_proj;
sys.krylov.Rpar_proj = U_0_sol;
end
% END OF FUNCTION AUX_PROPAGATION ----------------------

function [hom_tp,hom_tie] = aux_propagate_HomSol(c_sys,g_sys,timeStep,time,R0,taylorTerms,krylovError,dim_proj)
% propagate the homogeneous solution as well as the time-interval error
% stemming from the homogeneous solution

c = sparse(center(R0));
G = sparse(generators(R0));

% check if center is zero
c_norm = norm(c);

if isempty(c_sys)
    c_new = zeros(dim_proj,1);
    hom_tie = zeros(dim_proj,1);
else
    % Compute new center
    eAtk = readFieldForTimeStep(c_sys.taylor,'eAdt',time);
    if isempty(eAtk)
        eAtk = expm(c_sys.A*time);
        insertFieldTimeStep(c_sys.taylor,'eAdt',eAtk,time);
    end
    eAdt = readFieldForTimeStep(c_sys.taylor,'eAdt',timeStep);
    if isempty(eAdt)
        eAdt = expm(c_sys.A*timeStep);
        insertFieldTimeStep(c_sys.taylor,'eAdt',eAdt,timeStep);
    end
    c_expMatrix = eAdt * eAtk;
    insertFieldTimeStep(c_sys.taylor,'eAdt',c_expMatrix,time+timeStep);
    % c_expMatrix = expm(c_sys.A*time);
    c_new = c_norm*c_sys.B'*c_expMatrix(:,1);

     % Compute new TIE
    [~,F_c] = taylorMatrices(c_sys,timeStep,taylorTerms);
    hom_tie = c_norm*c_sys.B'*c_expMatrix*F_c(:,1);
end

% preallocation
nrOfGens = length(G(1,:));

% obtain generators using the Arnoldi iteration
G_new = zeros(dim_proj,size(G,2));

if nrOfGens > 0
    for iGen = 1:nrOfGens
    % parfor iGen = 1:nrOfGens
       g_norm = norm(G(:,iGen));
        if g_norm == 0
            G_new(:,iGen) = G(:,iGen);
        elseif g_norm <= 1e-8
            disp("HI");
        else
            %Compute new generator
            eAtk = readFieldForTimeStep(g_sys{iGen}.taylor,'eAdt',time);
            if isempty(eAtk)
                eAtk = expm(g_sys{iGen}.A*time);
                insertFieldTimeStep(g_sys{iGen}.taylor,'eAdt',eAtk,time);
            end
            eAdt = readFieldForTimeStep(g_sys{iGen}.taylor,'eAdt',timeStep);
            if isempty(eAdt)
                eAdt = expm(g_sys{iGen}.A*timeStep);
                insertFieldTimeStep(g_sys{iGen}.taylor,'eAdt',eAdt,timeStep);
            end
            expMatrix = eAdt * eAtk;
            insertFieldTimeStep(g_sys{iGen}.taylor,'eAdt',expMatrix,time+timeStep);
            % expMatrix = expm(g_sys{iGen}.A*time);
            G_new(:,iGen) = g_norm*g_sys{iGen}.B'*expMatrix(:,1);
        end
        % compute new TIE for this generator
        [~,F_g] = taylorMatrices(g_sys{iGen},timeStep,taylorTerms);
        hom_tie = hom_tie + g_norm*g_sys{iGen}.B'*expMatrix*F_g(:,1);        
    end
else
    G_new = []; %no generators 
end


if krylovError > 2*eps
    % Krylov error computation
    % +1 due to center
    error = krylovError * (size(R0.G,2) + 1);

    % initial-state-solution zonotope
    hom_tp = zonotope([c_new,G_new,error*eye(length(c_new))]);
    %----------------------------------------------------------------------
else
    hom_tp = zonotope([c_new,G_new]);
end


end
% END OF FUNCTION AUX_PROPAGATE_HOMSOL ------------------------------------


% ------------------------------ END OF CODE ------------------------------
