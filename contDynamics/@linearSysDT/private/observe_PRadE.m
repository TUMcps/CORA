function [R,tcomp] = observe_PRadE(linsysDT,params,options)
% observe_PRadE - computes the guaranteed state estimation approach
%    from [1].
%
% Syntax:
%    [R,tcomp] = observe_PRadE(linsysDT,params,options)
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
%    [1] Ye Wang, Meng Zhou, Vicenc Puig, Gabriela Cembrano, and
%        Zhenhua Wang. Zonotopic fault detection observer with H −
%        performance. In Proc. of the 36th IEEE Chinese Control
%        Conference, pages 7230–7235, 2017.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       18-September-2020
% Last update:   05-January-2021
%                25-February-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% obtain offline gains
options.L = observe_gain_PRadE(linsysDT,params,options);

% apply set-membership approach
tic
R = observe_intersectionFree(linsysDT,params,options);
tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
