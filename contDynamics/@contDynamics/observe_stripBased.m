function R = observe_stripBased(sys,params,options)
% observe_stripBased - computes the guaranteed state estimation approach
%    according to the set membership approach, see [1]
%
% Syntax:
%    R = observe_stripBased(sys,params,options)
%
% Inputs:
%    sys - continuous system object
%    params - model parameters
%    options - options for the computation of reachable sets
%
% Outputs:
%    R - estimated set of time points
%
% Reference:
%    [1] M. Althoff and J. J. Rath. Comparison of Set-Based Techniques 
%        for Guaranteed State Estimation of Linear Disturbed Systems. 
%        Automatica, 130, article no. 109662, 2021.
%    [2] M. Althoff. Guaranteed state estimation in CORA 2021. In Proc. 
%        of the 8th International Workshop on Applied Verification for 
%        Continuous and Hybrid Systems, 2021
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       08-September-2020
% Last update:   25-February-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%time period
tVec = params.tStart:options.timeStep:params.tFinal-options.timeStep;

% initialize parameter for the output equation
R = cell(length(tVec),1);

% width of strips
options.sigma = supremum(abs(interval(params.V)));

% Intersection
y = params.y(:,1);
% choose intersection procedure
if options.intersectionType == 1
    Rnext.tp = observe_intersectionMethod_I(sys,params.R0,y,options);
elseif options.intersectionType == 2
    Rnext.tp = observe_intersectionMethod_II(sys,params.R0,y,options);
end
% store first reachable set
R{1} = Rnext.tp;

% loop over all time steps
for k = 1:length(tVec)-1
    
    % Prediction
    Rnext = post(sys,Rnext,params.uTransVec(:,k),options);
    % add disturbance
    Rnext.tp = Rnext.tp + params.W;
    
    % Intersection
    y = params.y(:,k+1);
    % choose intersection procedure
    if options.intersectionType == 1
        Rnext.tp = observe_intersectionMethod_I(sys,Rnext.tp,y,options);
    elseif options.intersectionType == 2
        Rnext.tp = observe_intersectionMethod_II(sys,Rnext.tp,y,options);
    end
    
    % Order reduction
    Rnext.tp = reduce(Rnext.tp,options.reductionTechnique,options.zonotopeOrder);

    % Store result
    R{k+1} = Rnext.tp;
end

% ------------------------------ END OF CODE ------------------------------
