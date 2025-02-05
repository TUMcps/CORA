function res = test_interval_subsref
% test_interval_subsref - unit test function of subscripted reference
%
% Syntax:
%    res = test_interval_subsref
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

% Authors:       Dmitry Grebenyuk
% Written:       19-January-2016
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;

a = interval([-5.0; -4.0; -3; 0; 0; 5], [-2; 0.0; 2.0; 0; 5; 8]);
c = a;

assert( abs( infimum(c(1,1)) + 5.0 ) <= tol && abs( supremum(c(1,1)) + 2.0 ) <= tol )
assert( abs( infimum(c(2,1)) + 4.0 ) <= tol && abs( supremum(c(2,1)) + 0.0 ) <= tol )
assert( abs( infimum(c(3,1)) + 3.0 ) <= tol && abs( supremum(c(3,1)) - 2.0 ) <= tol )
assert( abs( infimum(c(4,1)) - 0.0 ) <= tol && abs( supremum(c(4,1)) - 0.0 ) <= tol )
assert( abs( infimum(c(5,1)) - 0.0 ) <= tol && abs( supremum(c(5,1)) - 5.0 ) <= tol )
assert( abs( infimum(c(6,1)) - 5.0 ) <= tol && abs( supremum(c(6,1)) - 8.0 ) <= tol )

% n-d arrays
inf = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
sup = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(inf,sup);
I_subs = I(:,:,1,:);
assert(isequal(size(I_subs),[2,2,1,3]));
I_subs = I(:,:,1); % index single page
assert(isequal(size(I_subs),[2,2]));

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
