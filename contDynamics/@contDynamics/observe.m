function [R,tcomp] = observe(sys,params,options)
% observe - computes the set of possible states of a set-based observer 
%
% Syntax:
%    R = observe(sys,params,options)
%
% Inputs:
%    sys - linearSysDT or nonlinearSysDT object
%    params - model parameters
%    options - options for bounding the set of states
%
% Outputs:
%    R - set of possible states of the observer
%    tcomp - computation time
%
%    [1] M. Althoff and J. J. Rath. Comparison of Set-Based Techniques 
%        for Guaranteed State Estimation of Linear Disturbed Systems. 
%        Automatica, 130, article no. 109662, 2021.
%    [2] M. Althoff. Guaranteed state estimation in CORA 2021. In Proc. 
%        of the 8th International Workshop on Applied Verification for 
%        Continuous and Hybrid Systems, 2021
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       14-June-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% currently only implemented for linearSysDT and nonlinearSysDT
if ~(isa(sys,'linearSysDT') || isa(sys,'nonlinearSysDT'))
    throw(CORAerror('CORA:noExactAlg',sys));
end

% options preprocessing
options = validateOptions(sys,mfilename,params,options);

% create vector of time points
tVec = (options.tStart:options.timeStep:options.tFinal-options.timeStep)';

% compute symbolic derivatives for nonlinear systems
if isa(sys,'nonlinearSysDT')
    derivatives(sys,options);
    % check output function of nonlinear system
    if ~all(sys.out_isLinear)
        throw(CORAerror('CORA:notSupported',...
            'Set-based observers only support linear output equations.'));
    end
end

% execute observer 
[R,tcomp] = executeObserver(sys,options);

% create object of class reachSet
timePoint.set = R;
timePoint.time = num2cell(tVec);
R = reachSet(timePoint);

% ------------------------------ END OF CODE ------------------------------
