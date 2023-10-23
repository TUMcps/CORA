function [R,tcomp] = observe_NomG(obj,options)
% observe_NomG - computes the guaranteed state estimation approach
% from [1].
%
% Syntax:
%    [R,tcomp] = observe_NomG(obj,options)
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
%    [1] Ye Wang, Vicenç Puig, and Gabriela Cembrano. Set-
%        membership approach and Kalman observer based on
%        zonotopes for discrete-time descriptor systems. Automatica,
%        93:435-443, 2018.
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

options.L = observe_gain_NomG(obj, options);

% apply set-membership approach
tic
R = observe_intersectionFree(obj, options);
tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
