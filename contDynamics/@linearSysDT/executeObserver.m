function [R,tcomp] = executeObserver(linsysDT,params,options)
% executeObserver - calls the appropriate observer
%
% Syntax:
%    [R,tcomp] = executeObserver(linsysDT,params,options)
%
% Inputs:
%    linsysDT - linearSysDT object
%    params - model parameters for guaranteed state estimation
%    options - options for the guaranteed state estimation
%
% Outputs:
%    R - observed set of points in time
%    tcomp - computation time
%
% References:
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
switch options.alg
    case 'VolMin-A' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_volMinA(linsysDT, params, options);
    case 'VolMin-B' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_volMinB(linsysDT, params, options);
    case 'FRad-A' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_FRadA(linsysDT, params, options);
    case 'FRad-B' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_FRadB(linsysDT, params, options);
    case 'PRad-A' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_PRadA(linsysDT, params, options);
    case 'PRad-B' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_PRadB(linsysDT, params, options);
    case 'PRad-C' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_PRadC(linsysDT, params, options);
    case 'FRad-C' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_FRadC(linsysDT, params, options);
    case 'PRad-D' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_PRadD(linsysDT, params, options);
    case 'PRad-E' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_PRadE(linsysDT, params, options);
    case 'Nom-G' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_NomG(linsysDT, params, options);
    case 'Hinf-G' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_HinfG(linsysDT, params, options);
    case 'ESO-A' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_ESO_A(linsysDT, params, options);
    case 'ESO-B' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_ESO_B(linsysDT, params, options);
    case 'ESO-C' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_ESO_C(linsysDT, params, options);    
    case 'ESO-D' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_ESO_D(linsysDT, params, options);
    case 'CZN-A' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_CZN_A(linsysDT, params, options);  
    case 'CZN-B' % set-based observer, see [1], [2]
        [R,tcomp] = priv_observe_CZN_B(linsysDT, params, options);  
    case 'ROPO' % set-based observer, see [3]
        [R,tcomp] = priv_observe_ROPO(linsysDT, params, options);
    case 'Rauch-Tung-Striebel' % smoother (not set-based), see [4]
        % here, the reachable set only consists of fixed values
        [R,tcomp] = priv_observe_RauchTungStriebel(linsysDT, params, options);
end

% ------------------------------ END OF CODE ------------------------------
