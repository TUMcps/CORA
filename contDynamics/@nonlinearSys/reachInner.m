function [Rin,Rout] = reachInner(nlnsys,params,options)
% reachInner - computes an inner-approximation of the reachable set
%
% Syntax:
%    Rin = reachInner(nlnsys,params,options)
%    [Rin,Rout] = reachInner(nlnsys,params,options)
%
% Inputs:
%    nlnsys - nonlinearSys object
%    params - parameter defining the reachability problem
%    options - struct containing the algorithm settings, crucially,
%              options.algInner = 'scale' -> algorithm in [1]
%              options.algInner = 'proj' -> algorithm in [2]
%              options.algInner = 'parallelo' -> algorithm in [3]
%              options.algInner = 'minkdiff' -> algorithm in [4]
%
% Outputs:
%    Rin - object of class reachSet storing the inner-approximation of the 
%          reachable set
%    Rout - object of class reachSet storing the outer-approximation of the
%           reachable set
%
% References:
%    [1] N. Kochdumper and M. Althoff. "Computing Non-Convex Inner-
%        Approximations of Reachable Sets for Nonlinear Continuous Systems"
%        CDC 2020
%    [2] E. Goubault and S. Putot. "Forward Inner-Approximated Reachability
%        of Non-Linear Continuous Systems", HSCC 2017 
%    [3] E. Goubault and S. Putot. "Robust Under-Approximations and 
%        Application to Reachability of Non-Linear Control Systems With 
%        Disturbances", Control System Letters 2021
%    [4] M. Wetzlinger, A. Kulmburg, and M. Althoff. "Inner approximations
%        of reachable sets for nonlinear systems using the Minkowski
%        difference". IEEE Control Systems Letters, 2024.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/reachInner

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       14-August-2020
% Last update:   17-December-2023 (MW, add MinkDiff algorithm)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% note: options preprocessing in subsequently called reachInner* functions

% compute inner-approximation with the selected algorithm
switch options.algInner

    % compute inner-approximation with the algorithm in [1]
    case 'scale'
        if nargout > 1
            [Rin,Rout] = priv_reachInnerScaling(nlnsys,params,options);
        else
            Rin = priv_reachInnerScaling(nlnsys,params,options);
        end

    % compute inner-approximation with the algorithm in [2]
    case 'proj'
        [Rin,Rout] = priv_reachInnerProjection(nlnsys,params,options);
        
    % compute inner-approximation with the algorithm in [3]    
    case 'parallelo'
        [Rin,Rout] = priv_reachInnerParallelotope(nlnsys,params,options);

    % compute inner-approximation with the algorithm in [4]
    case 'minkdiff'
        Rin = priv_reachInnerMinkdiff(nlnsys,params,options);
        Rout = [];

    otherwise
        throw(CORAerror('CORA:wrongFieldValue',...
            'options.algInner',{'scale','proj','parallelo','minkdiff'}));
end

% ------------------------------ END OF CODE ------------------------------
