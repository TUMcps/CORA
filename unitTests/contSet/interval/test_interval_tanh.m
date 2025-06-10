function res = test_interval_tanh
% test_interval_tanh - unit test function of tanh for intervals,
%    overloaded 'tanh()' function for intervals
%
% Syntax:
%    res = test_interval_tanh
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
I_tanh = tanh(I);
I_true = interval(-0.761594155955765,0.964027580075817);
assert(isequal(I_tanh,I_true,tol));

% unbounded
I = interval(-inf,2);
I_tanh = tanh(I);
I_true = interval(-1,0.964027580075817);
assert(isequal(I_tanh,I_true,tol));

I = interval(-1,inf);
I_tanh = tanh(I);
I_true = interval(-0.761594155955765,1);
assert(isequal(I_tanh,I_true,tol));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
