function [R,tcomp] = observe_ESO_A(obj,options)
% observe_ESO_A - computes the guaranteed state estimation approach
% from [1] and [2]. From [1] we use the intersection with strips and from
% [2] we use the prediction.
%
%
% Syntax:
%    [R,tcomp] = observe_ESO_A(obj,options)
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
%    [1] S. Gollamudi, S. Nagaraj, S. Kapoor, and Y. F. Huang.
%        Set-membership state estimation with optimal bounding
%        ellipsoids. In Proc. of the International Symposium on
%        Information Theory and its Applications, pages 262–265,
%        1996.
%    [2] Yushuang Liu, Yan Zhao, and Falin Wu. Ellipsoidal state-
%        bounding-based set-membership estimation for linear system
%        with unknown-but-bounded disturbances.
%        IET Control Theory & Applications, 10(4):431–442, 2016.
%
% Example: 
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
tVec = options.tStart:options.timeStep:options.tFinal-options.timeStep;

% initialize parameter for the output equation
R = cell(length(tVec),1);

% width of strips; here: use norm, see p.2 of [1]; todo: check
phi = norm(options.V);

% precompute result for eq. (48) in [2]
auxRes = sqrt(trace(options.W.Q));

% Intersection
y = options.y(:,1);
Rnext.tp = options.R0;
sigma_sq = 0.3;
% intersection of zonotope with strip
[Rnext.tp, sigma_sq] = intersectStrip(Rnext.tp,obj.C,...
    phi, y,sigma_sq,'Gollamudi1996');


% store first reachable set
R{1} = Rnext.tp;

% loop over all time steps
for k = 1:length(tVec)-1
    
    
    %% Prediction, eq. (21) of [2]
    % homogeneous solution
    Rnext_hom = post(obj,Rnext,options.uTransVec(:,k),options);
    % compute parameter p of eq. (48) in [2]
    p = auxRes/(sqrt(trace(sigma_sq*Rnext_hom.tp.Q))+auxRes);
    % finalize prediction
    Q_next = (1-p)^-1*Rnext_hom.tp.Q + (sigma_sq*p)^-1*options.W.Q;
    Rnext.tp = ellipsoid(Q_next, Rnext_hom.tp.q);
    
    % Intersection
    y = options.y(:,k+1);
    [Rnext.tp, sigma_sq] = intersectStrip(Rnext.tp,obj.C,...
        phi, y,sigma_sq,'Gollamudi1996');

    % Store result
    R{k+1} = Rnext.tp;
end
tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
