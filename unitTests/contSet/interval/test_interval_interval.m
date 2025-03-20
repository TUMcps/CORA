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

tol = 1e-12;

% empty interval
I = interval(zeros(4,0));
assert(representsa(I,'emptySet'))
assert(all(size(I.inf) == [4,0]) && all(size(I.sup) == [4,0]));

I = interval(zeros(4,0),zeros(4,0));
assert(representsa(I,'emptySet'))
assert(all(size(I.inf) == [4,0]) && all(size(I.sup) == [4,0]));

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
% standard case
I = interval(lb,ub);
assert(all(withinTol(I.inf,lb,tol)) && all(withinTol(I.sup,ub,tol)));

% point case
I = interval(lb);
assert(all(withinTol(I.inf,lb,tol)) && all(withinTol(I.sup,lb,tol)));

% point matrix case
I = interval(lb_mat);
assert(all(withinTol(I.inf,lb_mat,tol),"all") && all(withinTol(I.sup,lb_mat,tol),"all"));

% matrix case
I = interval(lb_mat,ub_mat);
assert(all(withinTol(I.inf,lb_mat,tol),"all") && all(withinTol(I.sup,ub_mat,tol),"all"));

% unbounded case
I = interval(lb_inf,ub_inf);
assert(all(withinTol(I.inf,lb_inf,tol)) && all(withinTol(I.sup,ub_inf,tol)));

% unbounded matrix case
I = interval(lb_mat_inf,ub_mat_inf);
assert(all(withinTol(I.inf,lb_mat_inf,tol),"all") && all(withinTol(I.sup,ub_mat_inf,tol),"all"));

% unbounded point case
I = interval(lb_inf,lb_inf);
assert(all(size(I.inf) == [0,2]) && all(size(I.sup) == [0,2]));

% unbounded point case 2
I = interval(ub_inf,ub_inf);
assert(all(size(I.inf) == [0,2]) && all(size(I.sup) == [0,2]));

% n-d arrays
lbnd = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ubnd = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lbnd,ubnd);
assert(all(withinTol(I.inf,lbnd),'all'))
assert(all(withinTol(I.sup,ubnd),'all'))

% combine results
res = true;


% wrong initializations
a_large = [10; 15];
b_small = [-20; -12];
a_plus1 = [-3; -5; -2; -8];
b_plus1 = [2; 6; 3; 9];

% empty instantiation
assertThrowsAs(@interval,'CORA:noInputInSetConstructor');

% lower limit larger than upper limit
assertThrowsAs(@interval,'CORA:wrongInputInConstructor',lb,b_small);
assertThrowsAs(@interval,'CORA:wrongInputInConstructor',a_large,ub);

% size of limits do not match
assertThrowsAs(@interval,'CORA:wrongInputInConstructor',a_plus1,ub);
assertThrowsAs(@interval,'CORA:wrongInputInConstructor',lb,b_plus1);
assertThrowsAs(@interval,'CORA:wrongInputInConstructor',lb_mat,ub);
assertThrowsAs(@interval,'CORA:wrongInputInConstructor',lb,ub_mat);


% too many input arguments
assertThrowsAs(@interval,'CORA:numInputArgsConstructor',lb,ub,ub);

% empty interval matrix via [-Inf,-Inf] or [Inf,Inf] entry
assertThrowsAs(@interval,'CORA:wrongInputInConstructor',lb_mat_inf);
assertThrowsAs(@interval,'CORA:wrongInputInConstructor',ub_mat_inf);

% ------------------------------ END OF CODE ------------------------------
