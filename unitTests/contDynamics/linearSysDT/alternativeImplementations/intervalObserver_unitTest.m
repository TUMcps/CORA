function estSet = intervalObserver_unitTest(obj,params,options)
% intervalObserver_unitTest - realizes the interval observer in [1].
%
% Syntax:
%    estSet = intervalObserver_unitTest(obj,params,options)
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
%    [1] W. Tang, Z. Wang, Y. Wang, T. Raissi, and Y. Shen.
%        Interval estimation methods for discrete-time linear time-
%        invariant systems. IEEE Transactions on Automatic Control,
%        64(11):4717-4724, 2019.

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


%% initialize
x = center(params.R0);
D = -options.L*params.V + params.W;
S_x = params.R0;
S_wv = zeros(length(x),1);

%% loop over all time steps
for k = 1:length(tVec)-1
    
    % error set
    E = box(S_x) + S_wv;
    
    % set of states
    S{k+1} = x + E;
    
    % update state
    x = obj.A*x + obj.B*params.u(:,k) + options.L*(params.y(:,k) - obj.C*x);

    % Update auxiliary sets
    S_x = (obj.A - options.L*obj.C)*S_x; 
    S_wv = S_wv + box(D);
    D = (obj.A - options.L*obj.C)*D;
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
