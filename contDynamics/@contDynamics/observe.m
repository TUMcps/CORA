function [R,tcomp] = observe(obj,params,options)
% observe - computes the set of possible states of a set-based observer 
%
% Syntax:  
%    R = observe(obj,params,options)
%
% Inputs:
%    obj - linearSysDT or nonlinearSysDT object
%    params - model parameters
%    options - options for bounding the set of states
%
% Outputs:
%    R - set of possible states of the observer
%    tcomp - computation time
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff
% Written:       14-Jun-2021
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% currently only implemented for linearSysDT and nonlinearSysDT
if ~(isa(obj,'linearSysDT') || isa(obj,'nonlinearSysDT'))
    throw(CORAerror('CORA:noExactAlg',obj));
end

% options preprocessing
options = validateOptions(obj,mfilename,params,options);

% create vector of time points
tVec = (options.tStart:options.timeStep:options.tFinal-options.timeStep)';

% compute symbolic derivatives for nonlinear systems
if isa(obj,'nonlinearSysDT')
    derivatives(obj,options);
    % check output function of nonlinear system
    if ~all(obj.out_isLinear)
        throw(CORAerror('CORA:notSupported',...
            'Set-based observers only support linear output equations.'));
    end
end

% execute observer 
[R,tcomp] = executeObserver(obj,options);

% create object of class reachSet
timePoint.set = R;
timePoint.time = num2cell(tVec);
R = reachSet(timePoint);

%------------- END OF CODE --------------