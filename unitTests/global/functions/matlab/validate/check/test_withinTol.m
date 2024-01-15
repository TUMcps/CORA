function res = test_withinTol
% test_withinTol - unit test function for checking equality of values
%    within tolerance
%
% Syntax:
%    res = test_withinTol
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
% Written:       30-April-2023
% Last update:   03-December-2023 (MW, add Inf cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% scalar values
a = 3;
b = 3.1;
res(end+1,1) = withinTol(a,a);
res(end+1,1) = ~withinTol(a,b);
res(end+1,1) = withinTol(a,b,b-a);

% scalar values (Inf)
a = -Inf;
b = Inf;
res(end+1,1) = withinTol(a,a);
res(end+1,1) = ~withinTol(a,b);

% scalar vs. vector
a = 0;
b = [0.001 0.002 -0.001 0 0.003];
res(end+1,1) = all(size(withinTol(a,b)) == size(b));
res(end+1,1) = all(withinTol(a,b,max(abs(b))-a));

% scalar vs. vector (Inf)
a = -Inf;
b = [-Inf;-Inf];
res(end+1,1) = all(size(withinTol(a,b)) == size(b));
res(end+1,1) = all(withinTol(a,b));

% vector vs. vector
a = [3 4 5];
b = [4 5 6];
res(end+1,1) = all(size(withinTol(a,b)) == size(a));
res(end+1,1) = ~any(withinTol(a,b));
res(end+1,1) = all(withinTol(a,b,1));

% scalar vs. matrix
a = 0;
B = [0 0.001 0.002; -0.001 -0.004 0.003];
res(end+1,1) = all(size(withinTol(a,B)) == size(B));
res(end+1,1) = all(all(withinTol(a,B,max(max(abs(B))))));

% scalar vs. matrix (Inf)
a = Inf;
B = Inf(3,2);
res(end+1,1) = all(size(withinTol(a,B)) == size(B));
res(end+1,1) = all(all(withinTol(a,B)));

% vector vs. matrix
a = [2 1];
B = [2 1.001; 1.999 0.999; 2.002 1.003];
res(end+1,1) = all(size(withinTol(a,B)) == size(B));
res(end+1,1) = all(all(withinTol(a,B,max(max(abs(B-a))))));

% vector vs. matrix (Inf)
a = [-Inf; Inf];
B = [-Inf -Inf; Inf Inf];
res(end+1,1) = all(size(withinTol(a,B)) == size(B));
res(end+1,1) = all(all(withinTol(a,B)));

% matrix vs. matrix
A = [0 1; 1 0; 2 0];
B = [-0.001 1; 1.001 0.002; 1.999 -0.002];
res(end+1,1) = all(size(withinTol(A,B)) == size(A));
res(end+1,1) = all(all(withinTol(A,B,max(max(abs(B-A))))));

% wrong syntaxes
try
    % too many input arguments
    withinTol(1,1,1,1);
    res = false;
end
try
    % vector of different lengths
    withinTol([1 0],[1 0 0]);
    res = false;
end
try
    % matrices of different size
    withinTol([1 0; 1 2],[1 0 0; 0 1 2]);
    res = false;
end
try
    % tolerance not a scalar
    withinTol([1 0],[1 0],[1 0]);
    res = false;
end
try
    % tolerance not a nonnegative value
    withinTol([1 0],[1 0],-1);
    res = false;
end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
