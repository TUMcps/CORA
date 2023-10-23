function [R,tcomp] = observe_PRadC(obj,options)
% observe_PRadC - computes the guaranteed state estimation approach
% from [1] and [2]. In [3], two versions of PRadC exist; we have removed
% PRadC-I since its performance is not as good as PRadC-II, which is now
% simply called PRad-C.
%
%
% Syntax:
%    [R,tcomp] = observe_PRadC(obj,options)
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
%    [2] Ye Wang, Teodoro Alamo, Vicenc Puig, and Gabriela
%        Cembrano. A distributed set-membership approach based on
%        zonotopes for interconnected systems. In Proc. of the IEEE
%        Conference on Decision and Control (CDC), pages 668–673, 2018.
%    [3] Althoff, M., Rath, J.~J.: Comparison of Guaranteed State 
%        Estimators for Linear Time-Invariant Systems , Automatica 130, 
%        2021, article no. 109662
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
% Last update:   04-January-2021
%                25-February-2021
%                01-July-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% obtain offline gains

OGain = observe_gain_PRadC(obj,options);

% set intersection procedure
options.intersectionType = 2;
options.intersectionTechnique = OGain; % gain directly provided

% apply set-membership approach
tic
R = observe_stripBased(obj,options);
tcomp = toc;

% ------------------------------ END OF CODE ------------------------------
