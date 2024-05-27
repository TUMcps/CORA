function res = testLong_interval_abs
% testLong_interval_abs - unit test function for the absolute value
%    computation of intervals
%
% Syntax:
%    res = testLong_interval_abs
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
% Written:       27-May-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% number of tests
nrTests = 1000;
res = true;
tol = eps;
nmax = 50;

for i=1:nrTests

    % random dimension(s): horizontal vector, vertical vector, matrix
    if mod(i,3) == 0
        n = [1, nmax];
    elseif mod(i,3) == 1
        n = [nmax, 1];
    else
        n = randi(nmax,1,2);
    end

    % lower bound larger than 0 -> abs has no effect
    lb = rand(n);
    ub = lb + 1;
    I = interval(lb,ub);
    I_abs = abs(I);

    if ~isequal(I,I_abs,tol)
        throw(CORAerror("CORA:testFailed"));
    end

    % upper bound smaller than 0 -> [-ub, -lb]
    ub = -rand(n);
    lb = ub - 1;
    I = interval(lb,ub);
    I_abs_true = interval(-ub, -lb);
    I_abs = abs(I);

    if ~isequal(I_abs_true,I_abs,tol)
        throw(CORAerror("CORA:testFailed"));
    end

    % lower and upper bound enclose 0, abs(lb) < abs(ub) -> [0, ub]
    lb = -rand(n);
    ub = -lb + 1;
    I = interval(lb,ub);
    I_abs_true = interval(zeros(n), ub);
    I_abs = abs(I);

    if ~isequal(I_abs_true,I_abs,tol)
        throw(CORAerror("CORA:testFailed"));
    end

    % lower and upper bound enclose 0, abs(lb) > abs(ub) -> [0, -lb]
    ub = rand(n);
    lb = -ub - 1;
    I = interval(lb,ub);
    I_abs_true = interval(zeros(n), -lb);
    I_abs = abs(I);
    if ~isequal(I_abs_true,I_abs,tol)
        throw(CORAerror("CORA:testFailed"));
    end

end

% ------------------------------ END OF CODE ------------------------------
