function res = test_interval_interval
% test_interval_interval - unit test function of interval constructor
%
% Syntax:
%    res = test_interval_interval
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
% Written:       27-July-2021
% Last update:   08-December-2023 (MW, add unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);
tol = 1e-12;

% empty interval
I = interval(zeros(4,0));
res(end+1,1) = representsa(I,'emptySet') && ...
    all(size(I.inf) == [4,0]) && all(size(I.sup) == [4,0]);

I = interval(zeros(4,0),zeros(4,0));
res(end+1,1) = representsa(I,'emptySet') && ...
    all(size(I.inf) == [4,0]) && all(size(I.sup) == [4,0]);

% random lower bound, random upper bound
lb = [-2; -3];
ub = [3; 5];
lb_inf = [-Inf, 2];
ub_inf = [3, Inf];
lb_mat = [-2 -1; 0 -3];
ub_mat = [3 2; 5 3];
lb_mat_inf = [-Inf -2; 2 4];
ub_mat_inf = [-5 -2; 3 Inf];

% admissible initializations
I = interval(lb,ub);
res(end+1,1) = all(withinTol(I.inf,lb,tol)) && all(withinTol(I.sup,ub,tol));

I = interval(lb);
res(end+1,1) = all(withinTol(I.inf,lb,tol)) && all(withinTol(I.sup,lb,tol));

I = interval(lb_mat);
res(end+1,1) = all(all(withinTol(I.inf,lb_mat,tol))) ...
    && all(all(withinTol(I.sup,lb_mat,tol)));

I = interval(lb_mat,ub_mat);
res(end+1,1) = all(all(withinTol(I.inf,lb_mat,tol))) ...
    && all(all(withinTol(I.sup,ub_mat,tol)));

I = interval(lb_inf,ub_inf);
res(end+1,1) = all(withinTol(I.inf,lb_inf,tol)) ...
    && all(withinTol(I.sup,ub_inf,tol));

I = interval(lb_mat_inf,ub_mat_inf);
res(end+1,1) = all(all(withinTol(I.inf,lb_mat_inf,tol))) ...
    && all(all(withinTol(I.sup,ub_mat_inf,tol)));

I = interval(lb_inf,lb_inf);
res(end+1,1) = all(size(I.inf) == [0,2]) && all(size(I.sup) == [0,2]);

I = interval(ub_inf,ub_inf);
res(end+1,1) = all(size(I.inf) == [0,2]) && all(size(I.sup) == [0,2]);

% combine results
res = all(res);


% wrong initializations
if CHECKS_ENABLED

a_large = [10; 15];
b_small = [-20; -12];
a_plus1 = [-3; -5; -2; -8];
b_plus1 = [2; 6; 3; 9];

% empty instantiation
try
    I = interval(); % <- should throw error here
    throw(CORAerror('CORA:testFailed'));
end

% lower limit larger than upper limit
try
    I = interval(lb,b_small); % <- should throw error here
    throw(CORAerror('CORA:testFailed'));
end
try
    I = interval(a_large,ub); % <- should throw error here
    throw(CORAerror('CORA:testFailed'));
end

% size of limits do not match
try
    I = interval(a_plus1,ub); % <- should throw error here
    throw(CORAerror('CORA:testFailed'));
end
try
    I = interval(lb,b_plus1); % <- should throw error here
    throw(CORAerror('CORA:testFailed'));
end
try
    I = interval(lb_mat,ub); % <- should throw error here
    throw(CORAerror('CORA:testFailed'));
end
try
    I = interval(lb,ub_mat); % <- should throw error here
    throw(CORAerror('CORA:testFailed'));
end

% too many input arguments
try
    I = interval(lb,ub,ub); % <- should throw error here
    throw(CORAerror('CORA:testFailed'));
end 

% empty interval matrix via [-Inf,-Inf] or [Inf,Inf] entry
try
    I = interval(lb_mat_inf); % <- should throw error here
    throw(CORAerror('CORA:testFailed'));
end

try
    I = interval(ub_mat_inf); % <- should throw error here
    throw(CORAerror('CORA:testFailed'));
end

end

% ------------------------------ END OF CODE ------------------------------
