function R = observe_intersectionFree(obj,options)
% observe_intersectionFree - computes the guaranteed state estimation 
% approach according to the intersection-free approach, see [1].
%
%
% Syntax:
%    R = observe_intersectionFree(obj,options)
%
% Inputs:
%    obj - discrete-time linear system object
%    options - options for the guaranteed state estimations
%
% Outputs:
%    R - observed set of points in time
%
% Reference:
%    [1] M. Althoff and J. J. Rath. Comparison of Set-Based Techniques 
%        for Guaranteed State Estimation of Linear Disturbed Systems, 
%        in preparation.
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       18-September-2020
% Last update:   05-January-2021
%                25-February-2021
%                26-February-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%%initialize computation

%time period
tVec = options.tStart:options.timeStep:options.tFinal-options.timeStep;
timeSteps = length(tVec);

% initialize parameter for the output equation
R = cell(length(tVec),1);

% store first reachable set
Rnext.tp = options.R0;
R{1} = Rnext.tp;


%% loop over all time steps
for k = 1:timeSteps-1
    
    % Prediction, eq. (11) in [1]
    Rnext.tp = (obj.A-options.L*obj.C)*Rnext.tp + obj.B*options.uTransVec(:,k) + ...
        options.L*options.y(:,k) + (-options.L*options.V) + options.W;
    
    % Order reduction
    Rnext.tp = reduce(Rnext.tp,options.reductionTechnique,options.zonotopeOrder);

    % Store result
    R{k+1} = Rnext.tp;
end


% ------------------------------ END OF CODE ------------------------------
