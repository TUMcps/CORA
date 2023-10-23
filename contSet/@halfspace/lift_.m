function hs = lift_(hs,N,dims)
% lift_ - lifts a halfspace object to a higher-dimensional space
%
% Syntax:
%    hs = lift_(hs,N,dims)
%
% Inputs:
%    hs - halfspace object
%    N - dimension of the higher dimensional space
%    dims - states of the high-dimensional space that correspond to the
%          states of the low-dimensional space
%
% Outputs:
%    hs - halfspace object in the higher-dimensional space
%
% Example: 
%    C = [2.1, 3.4, 5.2];
%    d = 1.7;
%    hs = halfspace(C,d);
%
%    hsHigh = lift(hs,10,[1,7,9])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/lift, polytope/lift, interval/lift

% Authors:       Niklas Kochdumper
% Written:       16-July-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize variables
C = zeros(N,1);
d = hs.d;

% insert parameters from the original halfspace object
C(dims) = hs.c;

% construct the resulting high dimensional halfspace object
hs = halfspace(C,d);

% ------------------------------ END OF CODE ------------------------------
