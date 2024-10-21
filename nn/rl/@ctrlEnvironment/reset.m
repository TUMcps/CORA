function [obj, initialObservation] = reset(obj)
% reset - reset ctrlEnvironment
%
% Syntax:
%   initialObservation = reset(obj)
%
% Inputs:
%   obj - ctrlEnvironment
% 
% Outputs:
%   obj - reset environment
%   initialObservation - initial observation
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ctrlEnvironment

% Authors:       Manuel Wendl
% Written:       03-November-2023
% Last update:   18-August-2024 (MW, new options structure)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

obj.stepNum = 1;

% Parse initial state for environment
if size(obj.options.rl.env.x0,1) ~= (obj.ctrlDynamics.nrOfOutputs-obj.ctrlDynamics.nrOfInputs)
    throw(CORAerror("CORA:dimensionMismatch",obj.options.rl.env.initialOps.x0,(obj.ctrlDynamics.nrOfOutputs-obj.ctrlDynamics.nrOfInputs)));
end

if strcmp(obj.options.rl.env.initialOps,'symmetric')
    r = binornd(1,0.5,size(obj.options.rl.env.x0));
    obj.state = obj.options.rl.env.x0.inf.*r + obj.options.rl.env.x0.sup.*~r;
elseif strcmp(obj.options.rl.env.initialOps,'uniform')
    obj.state = (obj.options.rl.env.x0.sup-obj.options.rl.env.x0.inf).*rand(size(obj.options.rl.env.x0))+obj.options.rl.env.x0.inf;
elseif strcmp(obj.options.rl.env.initialOps,'None')
    obj.state = 1/2*(obj.options.rl.env.x0.sup+obj.options.rl.env.x0.inf);
elseif strcmp(obj.options.rl.env.initialOps,'inf')
    obj.state = obj.options.rl.env.x0.inf;
elseif strcmp(obj.options.rl.env.initialOps,'sup')
    obj.state = obj.options.rl.env.x0.sup;
elseif strcmp(obj.options.rl.env.initialOps,'set') && strcmp(obj.options.rl.env.evalMode,'set')
    obj.state = obj.options.rl.env.x0;
else
    throw(CORAerror("CORA:notSupported",'This type of initial state seed is not supported'))
end

if strcmp(obj.options.rl.env.evalMode,'set')
    obj.state = zonotope(obj.state);
end

initialObservation = obj.state;

end

% ------------------------------ END OF CODE ------------------------------
