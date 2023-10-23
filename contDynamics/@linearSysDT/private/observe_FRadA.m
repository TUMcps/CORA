function [R,tcomp] = observe_FRadA(obj,options)
% observe_FRadA - computes the guaranteed state estimation approach
% from [1].
%
%
% Syntax:
%    [R,Rout] = observe_FRadA(obj,options)
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
%    [1] T. Alamo, J. M. Bravo, and E. F. Camacho. Guaranteed
%        state estimation by zonotopes. Automatica, 41(6):1035–1043,
%        2005.
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       12-September-2020
% Last update:   25-February-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set intersection procedure

options.intersectionType = 1;
options.intersectionTechnique = 'alamo-FRad';

% apply set-membership approach
tic
R = observe_stripBased(obj,options);
tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
