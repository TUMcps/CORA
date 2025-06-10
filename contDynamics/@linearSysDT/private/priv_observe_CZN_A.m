function [R,tcomp] = priv_observe_CZN_A(linsysDT,params,options)
% priv_observe_CZN_A - computes the guaranteed state estimation approach
%    from [1].
%
% Syntax:
%    [R,tcomp] = priv_observe_CZN_A(linsysDT,params,options)
%
% Inputs:
%    linsysDT - discrete-time linear system object
%    params - model parameters
%    options - options for the guaranteed state estimation
%
% Outputs:
%    R - observed set of points in time
%    tcomp - computation time
%
% Reference:
%    [1] J. K. Scott, D. M. Raimondo, G. R. Marseglia, and R. D.
%        Braatz. Constrained zonotopes: A new tool for set-based
%        estimation and fault detection. Automatica, 69:126–136, 2016.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       04-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

timerVal = tic;

% time period
tVec = params.tStart:options.timeStep:params.tFinal-options.timeStep;

% initialize parameter for the output equation
R = cell(length(tVec),1);

% Intersection; statement below eq. (33) of [1]
y = params.y(:,1);
Rnext.tp = conIntersect(params.R0, y + conZonotope(-1*params.V), linsysDT.C);

% store first reachable set
R{1} = Rnext.tp;

% loop over all time steps
for k = 1:length(tVec)-1
    
    % Prediction, part of eq. (33) of [1]
    Rnext.tp = linsysDT.A*Rnext.tp + linsysDT.B*params.uTransVec(:,k) + params.W;
    
    % Intersection, part of eq. (33) of [1]
    y = params.y(:,k+1);
    Rnext.tp = conIntersect(Rnext.tp, y + conZonotope(-1*params.V), linsysDT.C);
    
    % Order reduction
    Rnext.tp = reduce(Rnext.tp,options.reductionTechnique,options.zonotopeOrder);

    % Store result
    R{k+1} = Rnext.tp; % + 1e-6*zonotope([zeros(length(obj.C(1,:))), eye(length(obj.C(1,:)))]); % add small set for numerical stability
end

tcomp = toc(timerVal);

% ------------------------------ END OF CODE ------------------------------
