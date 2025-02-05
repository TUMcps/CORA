function [R,tcomp] = priv_observe_PRadA(linsysDT,params,options)
% priv_observe_PRadA - computes the guaranteed state estimation approach
%    from Sec. 4.1 in [1].
%
% Syntax:
%    [R,tcomp] = priv_observe_PRadA(linsysDT,params,options)
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
%    [1] V. T. H. Le, C. Stoica, T. Alamo, E. F. Camacho, and
%        D. Dumur. Zonotopic guaranteed state estimation for
%        uncertain systems. Automatica, 49(11):3418–3424, 2013.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       12-September-2020
% Last update:   02-January-2021
%                25-February-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% obtain offline gains
OGain = priv_observe_gain_PRadA(linsysDT,params,options);

% set intersection procedure
options.intersectionType = 1;
options.intersectionTechnique = OGain; % gain directly provided

% apply set-membership approach
tic
R = observe_stripBased(linsysDT,params,options);
tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
