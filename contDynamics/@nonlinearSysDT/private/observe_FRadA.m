function [R,tcomp] = observe_FRadA(nlnsysDT,params,options)
% observe_FRadA - computes the guaranteed state estimation approach
%    from [1]; the approach is extended here for nonlinear systems.
%
% Syntax:
%    [R,tcomp] = observe_FRadA(nlnsysDT,params,options)
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
%    [1] T. Alamo, J. M. Bravo, and E. F. Camacho. Guaranteed
%        state estimation by zonotopes. Automatica, 41(6):1035–1043, 2005.
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

% set intersection procedure
options.intersectionType = 1;
options.intersectionTechnique = 'alamo-FRad';

% apply set-membership approach
tic;
R = observe_stripBased(nlnsysDT,params,options);
tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
