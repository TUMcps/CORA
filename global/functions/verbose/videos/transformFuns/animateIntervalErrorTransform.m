function fun = animateIntervalErrorTransform(I)
% animateIntervalErrorTransform - returns a functions handle to animate 
%   adding an interval error
%
% Syntax:
%    res = animateIntervalErrorTransform(I)
%
% Inputs:
%    I - interval
%
% Outputs:
%    fun - function @(S,t) specifying transformed set S at time t \in [0,1]
%
% See also:
%    animateFromTo

% Authors:       Tobias Ladner
% Written:       16-December-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init function handle
fun = @(S,t) S + t*I;

end

% ------------------------------ END OF CODE ------------------------------
