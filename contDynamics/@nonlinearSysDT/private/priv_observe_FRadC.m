function [R,tcomp] = priv_observe_FRadC(nlnsysDT,params,options)
% priv_observe_FRadC - computes the guaranteed state estimation approach
%    from [1]; the approach is extended here for nonlinear systems.
%
% Syntax:
%    [R,tcomp] = priv_observe_FRadC(nlnsysDT,params,options)
%
% Inputs:
%    nlnsysDT - discrete-time nonlinear system object
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
% Written:       25-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set gain
options.observerType = 1; % combastel

% apply set-membership approach
tic;
R = priv_observe_intersectionFreeAdaptive(nlnsysDT,params,options);
tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
