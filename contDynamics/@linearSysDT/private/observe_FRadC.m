function [R,tcomp] = observe_FRadC(obj,options)
% observe_FRadC - computes the guaranteed state estimation approach
%    from [1].
%
% Syntax:  
%    [R,tcomp] = observe_FRadC(obj,options)
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
%    [1] Christophe Combastel. Zonotopes and kalman observers:
%        Gain optimality under distinct uncertainty paradigms and
%        robust convergence. Automatica, 55:265-273, 2015.
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff
% Written:       18-Sep-2020
% Last update:   25-Feb-2021
% Last revision: ---


%------------- BEGIN CODE --------------

% set gain
options.observerType = 1; % combastel

% apply set-membership approach
tic
R = observe_intersectionFreeAdaptive(obj,options);
tcomp = toc;

%------------- END OF CODE --------------
