function [R,tcomp] = observe_HinfG(obj,options)
% observe_HinfG - computes the guaranteed state estimation approach
% from [1].
%
%
% Syntax:
%    [R,tcomp] = observe_HinfG(obj,options)
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
%    [1] W. Tang, Z. Wang, Y. Wang, T. Raissi, and Y. Shen.
%        Interval estimation methods for discrete-time linear time-
%        invariant systems. IEEE Transactions on Automatic Control,
%        64(11):4717-4724, 2019.
%
% Example: 
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

options.L = observe_gain_HinfG(obj,options);

% apply set-membership approach
tic
R = observe_interval(obj,options);
tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
