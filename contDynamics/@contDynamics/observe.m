function [R, tcomp] = observe(obj,params,options)
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
    error("Function observe currently only implemented for discrete-time systems.");
end

% options preprocessing
options = validateOptions(obj,mfilename,params,options);

% create vector of time points
tVec = options.tStart:options.timeStep:options.tFinal - options.timeStep;
tVec = tVec';

% compute symbolic derivatives for nonlinear systems
if isa(obj,'nonlinearSysDT')
    derivatives(obj,options);
end

% execute observer 
[R,tcomp] = executeObserver(obj,options);

% create object of class reachSet
R = createReachSetObject(R,tVec,options);

end

% Auxiliary Functions -----------------------------------------------------

function R = createReachSetObject(R,tVec,options)
% create and object of class reachSet that stores the reachable set

    % time-point reachable set
    timePoint.set = R;
    timePoint.time = num2cell(tVec);
    
    if length(timePoint.time) ~= length(timePoint.set)
       timePoint.time = timePoint.time(1:length(timePoint.set)); 
    end
    
    % time interval reachable set (for plotting reasons, the time point set 
    % is also stored as a time interval set)
    timeInt.set = R;
    timeInt.time = cell(length(R),1);
    
    for i = 1:length(R)
        timeInt.time{i} = interval(tVec(i),tVec(i) + options.timeStep);
    end
    
    % construct object of class reachSet
    R = reachSet(timePoint, timeInt);
end


%------------- END OF CODE --------------