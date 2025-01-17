function fun = animateLinearTransform(varargin)
% animateLinearTransform - returns a functions handle to animate a linear
%   transformation
%
% Syntax:
%    res = animateLinearTransform(M,o)
%
% Inputs:
%    M - numeric matrix
%    o - numeric offset vector
%
% Outputs:
%    fun - function @(S,t) specifying transformed set S at time t \in [0,1]
%
% See also:
%    animateFromTo

% Authors:       Tobias Ladner
% Written:       11-December-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
[M,o] = setDefaultValues({1,0},varargin);

% idea: convex combination of identity In and matrix M
In = eye(size(M));

% init function handle
fun = @(S,t) ((1-t)*In + t*M) * S + t*o;

end

% ------------------------------ END OF CODE ------------------------------
