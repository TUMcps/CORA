function I = lift_(I,N,proj)
% lift_ - lifts an interval to a higher-dimensional space,
%
% Syntax:
%    S = lift_(S,N,proj)
%
% Inputs:
%    I - interval object
%    N - dimension of the higher-dimensional space
%    proj - states of the high-dimensional space that correspond to the
%          states of the low-dimensional interval object
%
% Outputs:
%    I - interval object in the higher-dimensional space
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/lift, interval/project, interval/projectHighDim

% Authors:       Tobias Ladner
% Written:       13-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init bounds
lb = -inf(N,1);
ub = inf(N,1);

% lift
lb(proj) = I.inf;
ub(proj) = I.sup;

% result
I = interval(lb, ub);

% ------------------------------ END OF CODE ------------------------------
