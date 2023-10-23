function estSet = setPropagationObserver_unitTest(obj,params,options)
% setPropagationObserver_unitTest - computes the guaranteed state 
% estimation approach according to the intersection-free approach
% when the gain changes in each iteration, see [1].
%
% Syntax:
%    estSet = setPropagationObserver_unitTest(obj,params,options)
%
% Inputs:
%    obj - linear discrete-time system object
%    params - model parameters
%    options - options for the computation
%
% Outputs:
%    estSet - set of estimated states
%
% Reference:
%    [1] C. Combastel. Zonotopes and Kalman observers:
%        Gain optimality under distinct uncertainty paradigms and
%        robust convergence. Automatica, 55:265-273, 2015.

% Authors:       ???
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%time period
tVec = params.tStart:options.timeStep:params.tFinal-options.timeStep;

% initialize set of estimated states
S = cell(length(tVec),1);

% store first set of estimated states
S{1} = params.R0;
Snext = S{1};


% loop over all time steps
for k = 1:length(tVec)-1

    % Prediction 
    Snext = (obj.A-options.L*obj.C)*Snext + obj.B*params.u(:,k) + ...
        options.L*params.y(:,k) + (-options.L*params.V) + params.W;
    
    % Order reduction
    Snext = reduce(Snext,options.reductionTechnique,options.zonotopeOrder);

    % Store result
    S{k+1} = Snext;
end

%% create object of class reachSet
% store sets
timeInt.set = S;

% store time intervals
timeInt.time = cell(length(S),1);
for i = 1:length(S)
    timeInt.time{i} = interval(tVec(i),tVec(i) + options.timeStep);
end

% construct object of class reachSet
estSet = reachSet([], timeInt);

% ------------------------------ END OF CODE ------------------------------
