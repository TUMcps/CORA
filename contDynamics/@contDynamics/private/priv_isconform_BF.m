function [res,R,failedDim] = priv_isconform_BF(sys,params,options)
% priv_isconform_BF - performs a conformance check using a white-box and test
%   data in a brute force manner. According to the definition of reachset 
%   conformance in [1], for each test case a reachability analysis is 
%   performed and it is checked whether the data points lie inside the 
%   reachable set. We only use the rechale outputs and not the reachable 
%   states. If no order reduction of the reachable set is performed, this 
%   method can be quite slow. The main purpose of this method is to test 
%   the correctness of more efficient methods.
%
% Syntax:
%    [res,R,failedDim] = priv_isconform_BF(sys,params,options)
%
% Inputs:
%    sys - contDynamics system
%    params - parameter defining the conformance problem
%    options - options for the conformance checking
%
% Outputs:
%    res - true if conformance was achieved, otherwise false
%    R - reachSet object (only time steps for which measurments exist)
%    failedDim - dimension in which enclosure failed; only applicable if
%                reachable set are enclosed by boxes.
%
% Reference:
%    [1] Hendrik Roehm, Jens Oehlerking, Matthias Woehrle, and Matthias 
%        Althoff. 2019. Model Conformance for Cyber-Physical Systems: 
%        A Survey. ACM Trans. Cyber-Phys. Syst. 3, 3, Article 30, 2019.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       07-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test suite
testSuite = params.testSuite;

% init partial results, reachable sets, and simulation results
resPartial = zeros(length(testSuite),1);
R = cell(length(testSuite),1);

% in case no initial set is provided
provideInitialState = true;

%% loop over test cases
for iCase = 1:length(testSuite)

    % compute time step and final time
    testCase = testSuite{iCase};
    timeSteps = size(testCase.y,1);
    options.timeStep = testCase.sampleTime/options.timeStepDivider;
    params.tFinal = testCase.sampleTime*(timeSteps-1);

    % input trajectory from test case
    for iSample = 1: size(testCase.u,3) % for linear systems
        params.u = testCase.u(:,:,iSample)';
        if isempty(params.u)
            % read default value in case none is given
            params.u = getDefaultValueParams('u',sys,params,options);
        end
        % consider time step divider
        n = size(params.u,1);
        tmp = [];
        for iStep = 1:options.timeStepDivider
            tmp((1:n) + (iStep-1)*n, :) = params.u;
        end
        params.u = reshape(tmp,n,[]);

        %% compute reachable set

        % rewrite params and options for reach
        [paramsReach,optionsReach] = aux_rewriteReach(params,options,provideInitialState,testCase);

        % compute reachable set
        R_full = reach(sys, paramsReach, optionsReach);

        % remove sets not part of the sampling time
        % init
        R_new.set = cell(timeSteps);
        R_new.time = cell(timeSteps);
        for iStep = 1:timeSteps
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
        R{iCase} = reachSet(R_new);

        %% enclosure check
        % initilaize failed dimensions
        failedDim.inf = [];
        failedDim.sup = [];
        % loop over all time steps
        for iStep = 1:timeSteps
            % check containment
            if ~contains(R{iCase}.timePoint.set{iStep}, testCase.y(iStep,:,iSample)')
                disp('Model is not reachset conformant');
                res = false;
                % return failed dimension for box enclosure
                if options.postProcessingOrder==1
                    I = interval(R{iCase}.timePoint.set{iStep});
                    % indices where infimum is breached
                    failedDim.inf = find(I.inf > testCase.y(iStep,:,iSample)');
                    % indices where supremum is breached
                    failedDim.sup = find(I.sup < testCase.y(iStep,:,iSample)');
                end
                return
            end
        end
    end
    
    % containment successfully shown
    resPartial(iCase) = 1;
end

% final conformance checking result
res = all(resPartial);

end


% Auxiliary functions -----------------------------------------------------

function [paramsReach,optionsReach] = aux_rewriteReach(params,options,provideInitialState,testCase)
    
    % copy params and options
    paramsReach = params;
    optionsReach = options;

    % update params ---
    % in case no initial set is provided
    if provideInitialState
        paramsReach.R0 = zonotope(testCase.initialState);
    end
    paramsReach = rmfield(paramsReach,'testSuite');


    % update options ---
    optionsReach.linAlg = options.reachAlg;
    optionsReach = rmiffield(optionsReach,'confAlg');
    optionsReach = rmiffield(optionsReach,'reachAlg');
    optionsReach = rmiffield(optionsReach,'timeStepDivider');
    optionsReach = rmiffield(optionsReach,'postProcessingOrder');
    optionsReach = rmiffield(optionsReach,'testSuite');
    optionsReach = rmiffield(optionsReach,'tFinal');
    optionsReach = rmiffield(optionsReach,'R0');
    optionsReach = rmiffield(optionsReach,'W');
    optionsReach = rmiffield(optionsReach,'V');
    optionsReach = rmiffield(optionsReach,'tStart');
    optionsReach = rmiffield(optionsReach,'U');
    optionsReach = rmiffield(optionsReach,'uTrans');
    optionsReach = rmiffield(optionsReach,'timeStep');
    optionsReach = rmiffield(optionsReach,'alg');
    optionsReach = rmiffield(optionsReach,'VALIDATED');

end

% ------------------------------ END OF CODE ------------------------------
