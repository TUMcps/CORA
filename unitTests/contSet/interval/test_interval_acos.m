function res = test_interval_acos
% test_interval_acos - unit test function of sine for intervals,
%    overloaded 'acos()' function for intervals
%
% Syntax:
%    res = test_interval_acos
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

% Authors:       Mark Wetzlinger
% Written:       08-August-2020
% Last update:   03-December-2023 (MW, add out of bounds cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;
res = true(0);

% Check special values:
% x     acos(x)
% -1    pi
% 0     pi/2
% 1     0
I = interval([-1,0],[1,0]);
I_acos = acos(I);
I_true = interval([0,pi/2],[pi,pi/2]);
res(end+1,1) = isequal(I_acos,I_true,tol);

% check out of bounds
I = interval(-2,0);
try
    I_acos = acos(I);
    res(end+1,1) = false;
end
I = interval(0.5,1.1);
try
    I_acos = acos(I);
    res(end+1,1) = false;
end


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
