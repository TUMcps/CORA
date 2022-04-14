function GSEstimation = evaluateSingleObserver(obj,params,options)
% evaluateSingleObserver - evaluates a specified observer
%
% Syntax:  
%    evaluateObservers()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Author:       Matthias Althoff
% Written:      08-Sep-2020
% Last update:  25-Feb-2021
% Last revision:---

%------------- BEGIN CODE --------------


% compute width of measurement strips
params.sigma = supremum(abs(interval(params.V)));

%% Run observer
[EstSet,tcomp] = observe(obj,params,options);

%% Performance Computation

% compute performance
perf = performanceObserver(EstSet);

% write all results to GSEstimation struct (guaranteed state estimation)
GSEstimation.Name = options.alg;
GSEstimation.ZOrder = options.zonotopeOrder;
GSEstimation.RedTech = options.reductionTechnique;
GSEstimation.EstStates = EstSet;
GSEstimation.TPlant = obj;
GSEstimation.InitialSet = params.R0;
GSEstimation.SampleStepSize = options.timeStep;
GSEstimation.Performance = perf;
GSEstimation.tIteration = tcomp/length(perf.volBox);


end