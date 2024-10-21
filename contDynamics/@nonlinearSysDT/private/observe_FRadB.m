function [R,tcomp] = observe_FRadB(nlnsysDT,params,options)
% observe_FRadB - computes the guaranteed state estimation approach
%    from [1]; the approach is extended here for nonlinear systems.
%
% Syntax:
%    [R,tcomp] = observe_FRadB(nlnsysDT,params,options)
%
% Inputs:
%    nlnsysDT - nonlinearSysDT object
%    params - model parameters
%    options - options for the guaranteed state estimation
%
% Outputs:
%    R - observed set of points in time
%    tcomp - computation time
%
% Reference:
%    [1] Ye Wang, Vicenç Puig, and Gabriela Cembrano. Set-membership
%        approach and Kalman observer based on zonotopes for discrete-time
%        descriptor systems. Automatica, 93:435-443, 2018.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       27-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set intersection procedure
options.intersectionType = 2;
options.intersectionTechnique.method = 'wang-FRad'; % type
options.intersectionTechnique.C = nlnsysDT.C; % measurement matrix
options.intersectionTechnique.E = generators(params.W); % disturbance matrix
options.intersectionTechnique.F = generators(params.V); % sensor noise matrix

% apply set-membership approach
tic;
R = observe_stripBased(nlnsysDT,params,options);
tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
