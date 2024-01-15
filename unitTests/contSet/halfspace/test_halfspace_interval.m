function res = test_halfspace_interval
% test_halfspace_interval - unit test function of interval conversion
%
% Syntax:
%    res = test_halfspace_interval
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       24-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty halfspace
n = 2;
hs = halfspace.empty(n);
I = interval(hs);
res(end+1,1) = isequal(I,interval.empty(n));


% 1D halfspaces
hs = halfspace(-3,4);
I = interval(hs);
res(end+1,1) = isequal(I,interval(-4/3,Inf));
hs = halfspace(2,5);
I = interval(hs);
res(end+1,1) = isequal(I,interval(-Inf,5/2));


% 2D halfspaces
hs = halfspace([2 0],1);
I = interval(hs);
res(end+1,1) = isequal(I,interval([-Inf;-Inf],[0.5;Inf]));
hs = halfspace([0 -3],2);
I = interval(hs);
res(end+1,1) = isequal(I,interval([-Inf;-2/3],[Inf;Inf]));
hs = halfspace([1 1],1);
I = interval(hs);
res(end+1,1) = isequal(I,interval([-Inf;-Inf],[Inf;Inf]));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
