function estSet = stripBasedObserver_unitTest(obj,params,options)
% stripBasedObserver_unitTest - computes the strip-based observer according to [1].
%
% Syntax:
%    estSet = stripBasedObserver_unitTest(obj,params,options)
%
% Inputs:
%    obj - discrete-time linear system object
%    params - model parameters
%    options - options for the guaranteed state estimation
%
% Outputs:
%    estSet - set of estimated states
%
% Reference:
%    [1] T. Alamo, J. M. Bravo, and E. F. Camacho. Guaranteed
%        state estimation by zonotopes. Automatica, 41(6):1035-1043,
%        2005.

% Authors:       ???
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%time period
tVec = params.tStart:options.timeStep:params.tFinal-options.timeStep;

% initialize parameter for the output equation
S = cell(length(tVec),1);

% width of strips
sigma = supremum(abs(interval(params.V)));

% store first set of estimated states
Snext = params.R0;

% Intersection of first measurement
y = params.y(:,1);
% loop over strips
for iStrip = 1:length(y)
    % intersection of zonotope with strip
    Snext = intersectStripAlamo_unitTest(Snext,obj.C(iStrip,:),sigma(iStrip),y(iStrip));
end

% Order reduction
Snext = reduce(Snext,options.reductionTechnique,options.zonotopeOrder);

% store first reachable set
S{1} = Snext;

% loop over all time steps
for k = 1:length(tVec)-1
    
    % Prediction
    Snext = obj.A*Snext + obj.B*params.u(:,k) + params.W;
    
    % Intersection
    y = params.y(:,k+1);
    % loop over strips
    for iStrip = 1:length(y)
        % intersection of zonotope with strip
        Snext = intersectStripAlamo_unitTest(Snext,obj.C(iStrip,:),sigma(iStrip),y(iStrip));
    end
    
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
