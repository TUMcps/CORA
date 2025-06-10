function [R,tcomp] = priv_observe_PRadB(linsysDT,params,options)
% priv_observe_PRadB - computes the guaranteed state estimation approach
%    from [1].
%
% Syntax:
%    [R,tcomp] = priv_observe_PRadB(linsysDT,params,options)
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
%        D. Dumur. Zonotope-based set-membership estimation for
%        multi-output uncertain systems. In Proc. of the IEEE
%        International Symposium on Intelligent Control (ISIC),
%        pages 212–217, 2013.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       18-September-2020
% Last update:   02-January-2021
%                25-February-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% obtain offline gains
OGain = priv_observe_gain_PRadB(linsysDT,params,options);

% set intersection procedure
options.intersectionType = 2;
options.intersectionTechnique = OGain; % gain directly provided

% apply set-membership approach
timerVal = tic;
R = observe_stripBased(linsysDT,params,options);
tcomp = toc(timerVal);

% ------------------------------ END OF CODE ------------------------------
