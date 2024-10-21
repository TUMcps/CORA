function [R,tcomp] = observe_ESO_B(linsysDT,params,options)
% observe_ESO_B - computes the guaranteed state estimation approach
%    from [1]
%
% Syntax:
%    [R,tcomp] = observe_ESO_B(linsysDT,params,options)
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
%    [1] Yushuang Liu, Yan Zhao, and Falin Wu. Ellipsoidal state-
%        bounding-based set-membership estimation for linear system
%        with unknown-but-bounded disturbances.
%        IET Control Theory & Applications, 10(4):431–442, 2016.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       03-March-2021
% Last update:   26-March-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tic

% time period
tVec = params.tStart:options.timeStep:params.tFinal-options.timeStep;

% initialize parameter for the output equation
R = cell(length(tVec),1);

% store auxiliary system information for set intersection
sys.V = params.V.Q;
sys.bar_V = chol(inv(params.V.Q)); % see  eq. (51) of [2]

% precompute result for eq. (48) in [1]
auxRes = sqrt(trace(params.W.Q));

% Intersection
y = params.y(:,1);
Rnext.tp = params.R0;
sigma_sq = 0.3;
% intersection of zonotope with strip
[Rnext.tp, sigma_sq] = intersectStrip(Rnext.tp,linsysDT.C,...
    sys, y, sigma_sq,'Liu2016');


% store first reachable set
R{1} = Rnext.tp;

% loop over all time steps
for k = 1:length(tVec)-1
    
    %% Prediction, eq. (21) of [1]
    % homogeneous solution
    Rnext_hom = post(linsysDT,Rnext,params.uTransVec(:,k),options);
    % compute parameter p of eq. (48) in [1]
    p = auxRes/(sqrt(trace(sigma_sq*Rnext_hom.tp.Q))+auxRes);
    % finalize prediction
    Q_next = (1-p)^-1*Rnext_hom.tp.Q + (sigma_sq*p)^-1*params.W.Q;
    Rnext.tp = ellipsoid(Q_next, Rnext_hom.tp.q, 1e-3);
    % Intersection
    y = params.y(:,k+1);
    [Rnext.tp, sigma_sq] = intersectStrip(Rnext.tp,linsysDT.C,...
        sys, y,sigma_sq,'Liu2016');

    % Store result
    R{k+1} = Rnext.tp;
end
tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
