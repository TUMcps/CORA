function [R,tcomp] = observe_FRadC(linsysDT,params,options)
% observe_FRadC - computes the guaranteed state estimation approach
%    from [1].
%
% Syntax:
%    [R,tcomp] = observe_FRadC(linsysDT,params,options)
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
%    [1] Christophe Combastel. Zonotopes and kalman observers:
%        Gain optimality under distinct uncertainty paradigms and
%        robust convergence. Automatica, 55:265-273, 2015.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       18-September-2020
% Last update:   25-February-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set gain
options.observerType = 1; % combastel

% apply set-membership approach
tic
R = observe_intersectionFreeAdaptive(linsysDT,params,options);
tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
