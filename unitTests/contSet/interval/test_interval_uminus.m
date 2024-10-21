function res = test_interval_uminus
% test_interval_uminus - unit test function of unary minus
%
% Syntax:
%    res = test_interval_uminus
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

% Authors:       Dmitry Grebenyuk, Mark Wetzlinger
% Written:       19-January-2016
% Last update:   04-December-2023 (MW, add empty and unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;

% empty
I = interval.empty(2);
I_uminus = -I;
assert(representsa(I_uminus,'emptySet'));

% bounded
I = interval([-5; -4; -3; 0; 0; 5], [-2; 0; 2; 0; 5; 8]);
I_uminus = -I;
I_true = interval([2;0;-2;0;-5;-8],[5;4;3;0;0;-5]);
assert(isequal(I_uminus,I_true,tol));

% unbounded
I = interval([-Inf;2],[1;Inf]);
I_uminus = -I;
I_true = interval([-1;-Inf],[Inf;-2]);
assert(isequal(I_uminus,I_true,tol));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
