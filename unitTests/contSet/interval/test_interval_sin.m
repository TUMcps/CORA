function res = test_interval_sin
% test_interval_sin - unit_test_function of sine for intervals,
%    overloaded 'sin()' function for intervals
%
% Syntax:
%    res = test_interval_sin
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

% Authors:       Dmitry Grebenyuk, Mark Wetzlinger
% Written:       13-January-2016
% Last update:   03-December-2023 (MW, add unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;
res = true;

% bounded
I = interval([-0; 0.0; 0.0; 0; 0; -pi/4; -2*pi], ...
             [pi/4; pi/2; pi; 2*pi; 4*pi; pi/4; 4*pi]);
I_sin = sin(I);
I_true = interval([0;0;0;-1;-1;-sqrt(2)/2;-1],[sqrt(2)/2;1;1;1;1;sqrt(2)/2;1]);
res(end+1,1) = isequal(I_sin,I_true,tol);

% unbounded
I = interval([-Inf;0;-Inf],[Inf;Inf;0]);
I_sin = sin(I);
I_true = interval(-ones(3,1),ones(3,1));
res(end+1,1) = isequal(I_sin,I_true,tol);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
