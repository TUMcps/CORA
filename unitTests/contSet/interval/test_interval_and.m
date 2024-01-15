function res = test_interval_and
% test_interval_and - unit test function of logical conjunction,
%    overloaded '&' operator for intervals
%
% Syntax:
%    res = test_interval_and
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
% See also: mtimes

% Authors:       Dmitry Grebenyuk, Mark Wetzlinger
% Written:       05-January-2016
% Last update:   23-April-2023 (MW, add empty set cases)
%                03-December-2023 (MW, add unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% bounded, bounded intersection
I1 = interval([-10; -2; 10], [-5; 8; 15]);
I2 = interval([-11; 2; 11], [-6; 9; 12]);
I = I1 & I2;
I_true = interval([-10;2;11],[-6;8;12]);
res(end+1,1) = isequal(I,I_true);

% bounded, degenerate intersection
I1 = interval([-2;3;1],[5;4;1]);
I2 = interval([5;3;0],[8;5;2]);
I = I1 & I2;
I_true = interval([5;3;1],[5;4;1]);
res(end+1,1) = isequal(I,I_true);

% bounded, empty intersection
I1 = interval([-5;-2],[2;4]);
I2 = interval([-7;6],[-3;8]);
I = I1 & I2;
res(end+1,1) = representsa(I,'emptySet');

% unbounded, bounded intersection
I1 = interval(-Inf,2);
I2 = interval(-2,Inf);
I = I1 & I2;
I_true = interval(-2,2);
res(end+1,1) = isequal(I,I_true);

% unbounded, unbounded intersection
I1 = interval(-Inf,2);
I2 = interval(-Inf,Inf);
I = I1 & I2;
I_true = interval(-Inf,2);
res(end+1,1) = isequal(I,I_true);

% unbounded, empty intersection
I1 = interval(-Inf,-2);
I2 = interval(2,Inf);
I = I1 & I2;
res(end+1,1) = representsa(I,'emptySet');

% unbounded, degenerate intersection
I1 = interval(-Inf,-1);
I2 = interval(-1,Inf);
I = I1 & I2;
I_true = interval(-1,-1);
res(end+1,1) = isequal(I,I_true);


% combine results 
res = all(res);

% ------------------------------ END OF CODE ------------------------------
