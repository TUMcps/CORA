function [R,tcomp] = observe_PRadD(obj,options)
% observe_PRadD - computes the guaranteed state estimation approach
% from [1].
%
%
% Syntax:
%    [R,tcomp] = observe_PRadD(obj,options)
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
%    [1] Ye Wang, Zhenhua Wang, Vicenc Puig, and Gabriela
%        Cembrano. Zonotopic set-membership state estimation for
%        discrete-time descriptor LPV systems. IEEE Transactions
%        on Automatic Control, 64(5):2092-2099, 2019.
%
% Example: 
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

OGain = observe_gain_PRadD(obj,options);

% set intersection procedure
options.intersectionType = 2;
options.intersectionTechnique = OGain; % gain directly provided

% apply set-membership approach
tic
R = observe_stripBased(obj,options);
tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
