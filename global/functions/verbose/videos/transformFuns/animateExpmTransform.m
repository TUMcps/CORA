function fun = animateExpmTransform(M,timeStep)
% animateExpmTransform - returns a functions handle to animate a 
%   matrix exponential transformation
%
% Syntax:
%    res = animateExpmTransform(M,timeStep)
%
% Inputs:
%    M - numeric matrix, e^M
%    timeStep - numeric, time step
%
% Outputs:
%    fun - function @(S,t) specifying transformed set S at time t \in [0,1]
%
% See also:
%    animateFromTo

% Authors:       Tobias Ladner
% Written:       17-December-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init function handle
fun = @(S,t) expm(M*(t*timeStep)) * S;

end

% ------------------------------ END OF CODE ------------------------------
