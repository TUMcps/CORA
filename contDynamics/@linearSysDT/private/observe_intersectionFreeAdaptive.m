function R = observe_intersectionFreeAdaptive(obj,options)
% observe_intersectionFreeAdaptive - computes the guaranteed state
% estimation approach according to the intersection-free approach
% when the gain changes in each iteration, see [1], [2].
%
% Syntax:
%    R = observe_intersectionFreeAdaptive(obj,options)
%
% Inputs:
%    obj - discrete-time linear system object
%    options - options for the guaranteed state estimation
%
% Outputs:
%    R - observed set of points in time
%
% Reference:
%    [1] Christophe Combastel. Zonotopes and kalman observers:
%        Gain optimality under distinct uncertainty paradigms and
%        robust convergence. Automatica, 55:265-273, 2015.
%    [2] M. Althoff and J. J. Rath. Comparison of Set-Based Techniques 
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
% Last update:   04-January-2021
%                25-February-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%time period

tVec = options.tStart:options.timeStep:options.tFinal-options.timeStep;

% initialize parameter for the output equation
R = cell(length(tVec),1);
% store first reachable set
Rnext.tp = options.R0;
R{1} = options.R0;

% F is chosen in [1] such that they are multiplied with unit
% uncertainties; thus, E and F can be seen as generators of zonotopes
% representing the disturbance and noise set
F = generators(options.V);

% loop over all time steps
for k = 1:length(tVec)-1
    
    % Compute observer gain
    if options.observerType == 1 % Combastel, FRad-C in [2]
        % obtain generators
        G = generators(Rnext.tp);
        G_comb = G*G';
%         L = obj.A * G_comb * obj.C' * inv(obj.C*G_comb*obj.C' + F*F');
        L = obj.A * G_comb * obj.C' / (obj.C*G_comb*obj.C' + F*F');
    else
        disp('this observer type is not yet implemented')
    end
    
    % Prediction, eq. (11) in [2]
    Rnext.tp = (obj.A-L*obj.C)*Rnext.tp + obj.B*options.uTransVec(:,k) + ...
        L*options.y(:,k) + (-L)*options.V + options.W;
    
    % Order reduction
    Rnext.tp = reduce(Rnext.tp,options.reductionTechnique,options.zonotopeOrder);
 
    % Store result
    R{k+1} = Rnext.tp;
end


% ------------------------------ END OF CODE ------------------------------
