function [R,tcomp] = observe_CZN_A(obj,options)
% observe_CZN_A - computes the guaranteed state estimation approach
% from [1].
%
%
% Syntax:
%    [R,tcomp] = observe_CZN_A(obj,options)
%
% Inputs:
%    obj - discrete-time linear system object
%    options - options for the guaranteed state estimation
%
% Outputs:
%    R - observed set of points in time
%    tcomp - computation time
%
% Reference:
%    [1] J. K. Scott, D. M. Raimondo, G. R. Marseglia, and R. D.
%        Braatz. Constrained zonotopes: A new tool for set-based
%        estimation and fault detection. Automatica, 69:126–136,
%        2016.
%
% Example: 
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

tic

% time period
tVec = options.tStart:options.timeStep:options.tFinal-options.timeStep;

% initialize parameter for the output equation
R = cell(length(tVec),1);

% Intersection; statement below eq. (33) of [1]
y = options.y(:,1);
Rnext.tp = conIntersect(options.R0, y + conZonotope(-1*options.V), obj.C);


% store first reachable set
R{1} = Rnext.tp;

% loop over all time steps
for k = 1:length(tVec)-1
    
    
    % Prediction, part of eq. (33) of [1]
    Rnext.tp = obj.A*Rnext.tp + obj.B*options.uTransVec(:,k) + options.W;
    
    % Intersection, part of eq. (33) of [1]
    y = options.y(:,k+1);
    Rnext.tp = conIntersect(Rnext.tp, y + conZonotope(-1*options.V), obj.C);
    
    % Order reduction
    Rnext.tp = reduce(Rnext.tp,options.reductionTechnique,options.zonotopeOrder);

    % Store result
    R{k+1} = Rnext.tp; % + 1e-6*zonotope([zeros(length(obj.C(1,:))), eye(length(obj.C(1,:)))]); % add small set for numerical stability
end
tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
