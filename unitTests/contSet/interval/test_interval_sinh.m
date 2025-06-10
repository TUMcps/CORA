function res = test_interval_sinh
% test_interval_sinh - unit test function of sinh for intervals,
%    overloaded 'sinh()' function for intervals
%
% Syntax:
%    res = test_interval_sinh
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

% Authors:       Tobias Ladner
% Written:       14-April-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;

% bounded
I = interval(-1,2);
I_sinh = sinh(I);
I_true = interval(-1.175201193643801,3.626860407847019);
assert(isequal(I_sinh,I_true,tol));

% unbounded
I = interval(-inf,2);
I_sinh = sinh(I);
I_true = interval(-Inf,3.626860407847019);
assert(isequal(I_sinh,I_true,tol));

I = interval(-1,inf);
I_sinh = sinh(I);
I_true = interval(-1.175201193643801,Inf);
assert(isequal(I_sinh,I_true,tol));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
