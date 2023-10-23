function [R,tcomp] = observe_FRadB(obj,options)
% observe_FRadB - computes the guaranteed state estimation approach
% from [1].
%
%
% Syntax:
%    [R,Rout] = observe_FRadB(obj,options)
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
% Written:       12-September-2020
% Last update:   04-January-2021
%                25-February-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set intersection procedure

options.intersectionType = 2;
options.intersectionTechnique.method = 'wang-FRad'; % type
options.intersectionTechnique.A = obj.A; % system matrix
options.intersectionTechnique.C = obj.C; % measurement matrix
options.intersectionTechnique.E = generators(options.W); % disturbance matrix
options.intersectionTechnique.F = generators(options.V); % sensor noise matrix

% apply set-membership approach
tic
R = observe_stripBased(obj,options);
tcomp = toc;

end

% ------------------------------ END OF CODE ------------------------------
