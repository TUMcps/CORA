function [Rin,Rout] = reachInner(sys,params,options)
% reachInner - compute an inner-approximation of the reachable set
%
% Syntax:
%    Rin = reachInner(sys,params,options)
%    [Rin,Rout] = reachInner(sys,params,options)
%
% Inputs:
%    sys - nonlinearSys object
%    params - parameter defining the reachability problem
%    options - struct containing the algorithm settings
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
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/reachInner

% Authors:       Niklas Kochdumper
% Written:       14-August-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% note: options preprocessing in reachInner* functions

% compute inner-approximation with the selected algorithm
switch options.algInner

    % compute inner-approximation with the algorithm in [1]
    case 'scale'
        if nargout > 1
            [Rin,Rout] = reachInnerScaling(sys,params,options);
        else
            Rin = reachInnerScaling(sys,params,options);
        end

    % compute inner-approximation with the algorithm in [2]
    case 'proj'
        [Rin,Rout] = reachInnerProjection(sys,params,options);
        
    % compute inner-approximation with the algorithm in [3]    
    case 'parallelo'
        [Rin,Rout] = reachInnerParallelotope(sys,params,options);

    otherwise
        throw(CORAerror('CORA:wrongFieldValue',...
            'options.algInner',{'scale','proj','parallelo'}));
end

% ------------------------------ END OF CODE ------------------------------
