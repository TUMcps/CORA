function fun = animateNoTransform()
% animateNoTransform - returns a functions handle to animate nothing.
%   This "transformation" usually speeds up the animation 
%   if only colors need to be adjusted.
%
% Syntax:
%    res = animateNoTransform()
%
% Inputs:
%    -
%
% Outputs:
%    fun - function @(S,t) specifying transformed set S at time t \in [0,1]
%
% See also:
%    animateFromTo

% Authors:       Tobias Ladner
% Written:       12-December-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init function handle
fun = @(S,t) S;

end

% ------------------------------ END OF CODE ------------------------------
