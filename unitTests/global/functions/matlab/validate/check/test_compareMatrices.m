function res = test_compareMatrices
% test_compareMatrices - unit test function for comparsion of matrices
%
% Syntax:
%    res = test_compareMatrices()
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
% Written:       22-November-2022
% Last update:   08-May-2023 (TL, ordered)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true, wait for failure
res = true;

% init matrices
A = [1 2 3 4; 5 6 7 8; 9 10 11 12];
B = [1 2 3 4; 5 6 7 8; 9 10 11 12];
B_order = B(:,[4 2 3 1]);
B_eps = B_order + 10*eps;
B_minus = B_order(:,1:3);
B_plus = [B_order, [1; 1; 1]];
C = [9 10 11 12; 5 6 7 8; 1 2 3 4];
D = [];

% completely equal matrices
if ~compareMatrices(A,B)
    res = false;
elseif ~compareMatrices(A,B,eps,'subset')
    res = false;
end

% compare with different order
if ~compareMatrices(A,B_order)
    res = false;
elseif ~compareMatrices(A,B_order,eps,'subset')
    res = false;
end

% check with tolerance
if compareMatrices(A,B_eps,eps)
    % tolerance too small
    res = false;
elseif ~compareMatrices(A,B_eps,100*eps)
    % tolerance large enough
    res = false;
end

% check different sizes and subset
if compareMatrices(A,B_plus)
    res = false;
elseif compareMatrices(A,B_minus)
    res = false;
end

% check subset
if ~compareMatrices(B_minus,B,eps,'subset')
    res = false;
elseif ~compareMatrices(A,B_plus,eps,'subset')
    res = false;
end

% completely different matrices
if compareMatrices(A,C)
    res = false;
end

% check ordered
if ~compareMatrices(A,B,eps,'equal',true)
    res = false;
elseif compareMatrices(A(:,1:end-1),B,eps,'equal',true)
    res = false;
elseif ~compareMatrices(A,B,eps,'subset',true)
    res = false;
elseif ~compareMatrices(A(:,1:end-1),B,eps,'subset',true)
    res = false;
elseif compareMatrices(A,B_order,eps,'equal',true)
    res = false;
elseif compareMatrices(A,B_order,eps,'subset',true)
    res = false;
end


% check wrong input arguments (should not reach the next line)
if CHECKS_ENABLED

try
    % matrices have to be nonempty
    compareMatrices(A,D);
    res = false;
end
try
    % matrices have to be nonempty
    compareMatrices(D,B);
    res = false;
end
try
    % eps has to be nonnegative
    compareMatrices(A,B,-1);
    res = false;
end
try
    % eps has to be a scalar
    compareMatrices(A,B,[1 1]);
    res = false;
end
try
    % flag has to be fourth input argument
    compareMatrices(A,B,'subset');
    res = false;
end
try
    % flag has to be 'equal' or 'subset'
    compareMatrices(A,B,eps,'exact');
    res = false;
end

end

% ------------------------------ END OF CODE ------------------------------
