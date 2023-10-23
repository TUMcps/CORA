function [R,tcomp] = executeObserver(obj,options)
% executeObserver - calls the appropriate observer
%
% Syntax:
%    [R,tcomp] = executeObserver(obj,options)
%
% Inputs:
%    obj - discrete-time linear system object
%    options - options for the guaranteed state estimation
%
% Outputs:
%    R - observed set of points in time
%    tcomp - computation time
%
%    [1] M. Althoff and J. J. Rath. Comparison of Set-Based Techniques 
%        for Guaranteed State Estimation of Linear Disturbed Systems. 
%        Automatica, 130, article no. 109662, 2021.
%    [2] M. Althoff. Guaranteed state estimation in CORA 2021. In Proc. 
%        of the 8th International Workshop on Applied Verification for 
%        Continuous and Hybrid Systems, 2021
%    [3] Vicino, A., & Zappa, G. (1996). Sequential approximation of 
%        feasible parameter sets for identification with set membership 
%        uncertainty. IEEE Transactions on Automatic Control, 41(6), 774-785.
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       20-March-2020
% Last update:   25-February-2021
%                14-June-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% decide which observer to execute by options.alg
if strcmp(options.alg,'VolMin-A') % set-based observer, see [1], [2]
    [R,tcomp] = observe_volMinA(obj, options);
elseif strcmp(options.alg,'VolMin-B') % set-based observer, see [1], [2]
    [R,tcomp] = observe_volMinB(obj, options);
elseif strcmp(options.alg,'FRad-A') % set-based observer, see [1], [2]
    [R,tcomp] = observe_FRadA(obj, options);
elseif strcmp(options.alg,'FRad-B') % set-based observer, see [1], [2]
    [R,tcomp] = observe_FRadB(obj, options);
elseif strcmp(options.alg,'PRad-A') % set-based observer, see [1], [2]
    [R,tcomp] = observe_PRadA(obj, options);
elseif strcmp(options.alg,'PRad-B') % set-based observer, see [1], [2]
    [R,tcomp] = observe_PRadB(obj, options);
elseif strcmp(options.alg,'PRad-C') % set-based observer, see [1], [2]
    [R,tcomp] = observe_PRadC(obj, options);
elseif strcmp(options.alg,'FRad-C') % set-based observer, see [1], [2]
    [R,tcomp] = observe_FRadC(obj, options);
elseif strcmp(options.alg,'PRad-D') % set-based observer, see [1], [2]
    [R,tcomp] = observe_PRadD(obj, options);
elseif strcmp(options.alg,'PRad-E') % set-based observer, see [1], [2]
    [R,tcomp] = observe_PRadE(obj, options);
elseif strcmp(options.alg,'Nom-G') % set-based observer, see [1], [2]
    [R,tcomp] = observe_NomG(obj, options);
elseif strcmp(options.alg,'Hinf-G') % set-based observer, see [1], [2]
    [R,tcomp] = observe_HinfG(obj, options);
elseif strcmp(options.alg,'ESO-A') % set-based observer, see [1], [2]
    [R,tcomp] = observe_ESO_A(obj, options);
elseif strcmp(options.alg,'ESO-B') % set-based observer, see [1], [2]
    [R,tcomp] = observe_ESO_B(obj, options);
elseif strcmp(options.alg,'ESO-C') % set-based observer, see [1], [2]
    [R,tcomp] = observe_ESO_C(obj, options);    
elseif strcmp(options.alg,'ESO-D') % set-based observer, see [1], [2]
    [R,tcomp] = observe_ESO_D(obj, options);
elseif strcmp(options.alg,'CZN-A') % set-based observer, see [1], [2]
    [R,tcomp] = observe_CZN_A(obj, options);  
elseif strcmp(options.alg,'CZN-B') % set-based observer, see [1], [2]
    [R,tcomp] = observe_CZN_B(obj, options);  
elseif strcmp(options.alg,'ROPO') % set-based observer, see [3]
    [R,tcomp] = observe_ROPO(obj,options);
elseif strcmp(options.alg,'Rauch-Tung-Striebel') % smoother (not set-based), see [4]
    [R,tcomp] = observe_RauchTungStriebel(obj,options); % here, the reachable set only consists of fixed values
end

% ------------------------------ END OF CODE ------------------------------
