function res = test_polytope_interval
% test_polytope_interval - unit test function of conversion to intervals
%
% Syntax:
%    res = test_polytope_interval
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
% Written:       29-November-2022
% Last update:   14-December-2022 (MW, add unbounded cases)
%                27-July-2023 (MW, 1D cases) 
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D, unbounded
P = polytope(1,1);
I = interval(P);
I_true = interval(-Inf,1);
res(end+1,1) = isequal(I,I_true);

% 1D, bounded
P = polytope([3;-2],[1;0]);
I = interval(P);
I_true = interval(0,1/3);
res(end+1,1) = isequal(I,I_true);

% 1D, single point
P = polytope([],[],5,3);
I = interval(P);
I_true = interval(3/5);
res(end+1,1) = isequal(I,I_true);


% 2D, diamond in unit square
P = polytope([1 1; -1 1; -1 -1; 1 -1],ones(4,1));
% convert to interval and compare to true result
I = interval(P);
I_true = interval([-1;-1],[1;1]);
res(end+1,1) = isequal(I,I_true);


% 3D, unit cube
n = 3;
P = polytope([eye(n); -eye(n)],ones(2*n,1));
% convert to interval and show exact conversion
I = interval(P);
res(end+1,1) = isequal(P,I);


% 4D, init unbounded polytope
A = [1 0 0 0; 0 1 0 0; 0 0 0 1; -1 0 0 0; 0 0 0 -1];
b = [1; 2; 4; 3; 1.5];
P = polytope(A,b);
% convert to interval and compare to true results
I = interval(P);
I_true = interval([-3; -Inf; -Inf; -1.5], [1; 2; Inf; 4]);
res(end+1,1) = isequal(I,I_true);


% combine results
res = all(res);


% ------------------------------ END OF CODE ------------------------------
