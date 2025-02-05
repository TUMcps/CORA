function [R,tcomp] = priv_observe_NomG(linsysDT,params,options)
% priv_observe_NomG - computes the guaranteed state estimation approach
%    from [1].
%
% Syntax:
%    [R,tcomp] = priv_observe_NomG(linsysDT,params,options)
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
%    [1] Ye Wang, Vicenç Puig, and Gabriela Cembrano. Set-
%        membership approach and Kalman observer based on
%        zonotopes for discrete-time descriptor systems. Automatica,
%        93:435-443, 2018.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       05-January-2021
% Last update:   25-February-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% obtain offline gains
options.L = priv_observe_gain_NomG(linsysDT, params, options);

% apply set-membership approach
tic
R = priv_observe_intersectionFree(linsysDT, params, options);
tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
