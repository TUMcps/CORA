function [res,R,traj_RRT,failedDim] = priv_isconform_RRT(sys,params,options)
% priv_isconform_RRT - performs a conformance check using a white-box model and
%   rapidly exploring random trees (RRTs), see [1].
%
% Syntax:
%    [res,R,traj,failedDim] = priv_isconform_RRT(sys,params,options)
%
% Inputs:
%    sys - contDynamics system
%    params - parameter defining the conformance problem
%    options - options for the conformance checking
%
% Outputs:
%    res - true if conformance was achieved, otherwise false
%    R - reachSet object (only time steps for which measurments exist)
%    traj - states of the rapidly exploring random tree
%    failedDim - dimension for which the conformance check failed (only for box enclosure of reachable sets)
%
% Reference:
%    [1] M. Althoff and J. M. Dolan. Reachability computation of low-order 
%        models for the safety verification of high-order road vehicle 
%        models. In Proc. of the American Control Conference, 
%        page 3559–3566, 2012.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       16-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test suite
testSuite = params.testSuite;

% init partial results, reachable sets, and simulation results
resPartial = zeros(length(testSuite),1);
R = cell(length(testSuite),1);
traj_RRT = [];

%% loop over trajectories
for m = 1:length(testSuite)

    % compute time step and final time
    traj_m = testSuite(m);
    timeSteps = size(traj_m.u,2);
    options.timeStep = traj_m.dt/options.timeStepDivider;
    params.tFinal = traj_m.dt*timeSteps;
    
    
    % input trajectory from test case
    params.u = traj_m.u;
    % consider time step divider
    dim = size(traj_m.u,1);
    tmp = [];
    for iStep = 1:options.timeStepDivider
        tmp((1:dim) + (iStep-1)*dim, :) = params.u;
    end
    params.u = reshape(tmp,dim,[]);

    %% compute reachable set

    % rewrite params and options for reach
    [paramsReach,optionsReach] = aux_rewriteReach(params,options);

    % compute reachable set
    R_full = reach(sys, paramsReach, optionsReach);
    % remove sets not part of the sampling time
    % init
    R_new.set = cell(timeSteps+1);
    R_new.time = cell(timeSteps+1);
    for iStep = 1:timeSteps+1
        % original index
        ind = (iStep-1)*options.timeStepDivider + 1;
        % saving using new index
        R_new.set{iStep} = R_full.timePoint.set{ind};
        R_new.time{iStep} = R_full.timePoint.time{ind};
        % post-processing order reduction (has to be performed equally to
        % verify the system)
        R_new.set{iStep} = reduce(R_new.set{iStep},options.reductionTechnique,options.postProcessingOrder);
    end
    % init reachset
    R{m} = reachSet(R_new);

    %% compute rapidly-exploring random tree
    % "hack" before validateOptions is fixed (remove this later)
    % rewrite to remove extended options.uTransVec for intermediate time
    % steps in reachability analysis
    options.uTransVec = traj_m.u;
    % check for pre-computed RRT
    if isfield(options,'preComputedRRT')
        load(options.preComputedRRT, 'traj_RRT')
    else
        % compute RRT
        options.type = 'rrt';
        options.R = R{m};
        traj_RRT = [trajRRT; simulateRandom(options.refModel, params, options)];
        % save result
        save trajRRT traj_RRT
    end

    %% enclosure check
    % initilaize failed dimensions
    failedDim.inf = [];
    failedDim.sup = [];
    % loop over all time steps
    for iStep = 1:timeSteps+1
        % loop over all samles
        for iSample = 1:options.points
            % project result form high-fidelity model
            xProj = options.convertToAbstractState(traj_RRT(m).x(:,iStep,iSample));
            % check containment
            if ~contains(R{m}.timePoint.set{iStep}, xProj)
                disp('Model is not reachset conformant');
                res = false;
                % return failed dimension for box enclosure
                if options.postProcessingOrder==1
                    I = interval(R{m}.timePoint.set{iStep});
                    % indices where infimum is breached
                    failedDim.inf = find(I.inf > xProj);
                    % indices where supremum is breached
                    failedDim.sup = find(I.sup < xProj);
                end
                return
            end
        end
    end
    
    % containment successfully shown
    resPartial(m) = 1;

end

% final conformance checking result
res = all(resPartial);

end


% Auxiliary functions -----------------------------------------------------

function [paramsReach,optionsReach] = aux_rewriteReach(params,options)
    
    % copy params and options
    paramsReach = params;
    optionsReach = options;

    % update params ---
    paramsReach = rmiffield(paramsReach,'testSuite');

    % update options ---
    optionsReach.alg = options.reachAlg;
    optionsReach = rmiffield(optionsReach,'timeStepDivider');
    optionsReach = rmiffield(optionsReach,'postProcessingOrder');
    optionsReach = rmiffield(optionsReach,'confAlg');
    optionsReach = rmiffield(optionsReach,'reachAlg');
    optionsReach = rmiffield(optionsReach,'U_increment');
    optionsReach = rmiffield(optionsReach,'points');
    optionsReach = rmiffield(optionsReach,'vertSamp');
    optionsReach = rmiffield(optionsReach,'stretchFac');
    optionsReach = rmiffield(optionsReach,'convertFromAbstractState');
    optionsReach = rmiffield(optionsReach,'convertToAbstractState');
    optionsReach = rmiffield(optionsReach,'preComputedRRT');
    optionsReach = rmiffield(optionsReach,'refModel');
    optionsReach = rmiffield(optionsReach,'uTransVec');
    % optionsReach = rmiffield(optionsReach,'intermediateOrder');

end

% ------------------------------ END OF CODE ------------------------------
