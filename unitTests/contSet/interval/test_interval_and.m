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

% bounded, bounded intersection
I1 = interval([-10; -2; 10], [-5; 8; 15]);
I2 = interval([-11; 2; 11], [-6; 9; 12]);
I = I1 & I2;
I_true = interval([-10;2;11],[-6;8;12]);
assert(isequal(I,I_true));

% bounded, degenerate intersection
I1 = interval([-2;3;1],[5;4;1]);
I2 = interval([5;3;0],[8;5;2]);
I = I1 & I2;
I_true = interval([5;3;1],[5;4;1]);
assert(isequal(I,I_true));

% bounded, empty intersection
I1 = interval([-5;-2],[2;4]);
I2 = interval([-7;6],[-3;8]);
I = I1 & I2;
assert(representsa(I,'emptySet'));

% unbounded, bounded intersection
I1 = interval(-Inf,2);
I2 = interval(-2,Inf);
I = I1 & I2;
I_true = interval(-2,2);
assert(isequal(I,I_true));

% unbounded, unbounded intersection
I1 = interval(-Inf,2);
I2 = interval(-Inf,Inf);
I = I1 & I2;
I_true = interval(-Inf,2);
assert(isequal(I,I_true));

% unbounded, empty intersection
I1 = interval(-Inf,-2);
I2 = interval(2,Inf);
I = I1 & I2;
assert(representsa(I,'emptySet'));

% unbounded, degenerate intersection
I1 = interval(-Inf,-1);
I2 = interval(-1,Inf);
I = I1 & I2;
I_true = interval(-1,-1);
assert(isequal(I,I_true));

% empty intervals
I = interval([1;2],[3;4]);
Iempty = interval.empty(2);
I_and = I & Iempty;
assert(isequal(Iempty,I_and));
I_and = Iempty & I;
assert(isequal(Iempty,I_and));

% n-d arrays
lb = [];
lb(:,:,1,1) = [1 2; 3 5];
lb(:,:,1,2) = [0 -1; -2 3];
lb(:,:,1,3) = [1 1; -1 0];
lb(:,:,2,1) = [-3 2; 0 1];
ub = [];
ub(:,:,1,1) = [1.5 4; 4 10];
ub(:,:,1,2) = [1 2; 0 4];
ub(:,:,1,3) = [2 3; -0.5 2];
ub(:,:,2,1) = [-1 3; 0 2];
I = interval(lb,ub);
Iabs = abs(I);
assert(representsa(and(I,Iabs),'emptySet'));

% combine results 
res = true;

% ------------------------------ END OF CODE ------------------------------
