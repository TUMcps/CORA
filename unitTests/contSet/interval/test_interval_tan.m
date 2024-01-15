function res = test_interval_tan
% test_interval_tan - unit test function of tangent function
%
% Syntax:
%    res = test_interval_tan
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

% Authors:       Daniel Althoff, Mark Wetzlinger
% Written:       03-November-2015
% Last update:   04-January-2016
%                04-December-2023 (MW, update syntax, add unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-9;
res = true(0);

% bounded, degenerate
I = interval(0);
I_tan = tan(I);
res(end+1,1) = isequal(I_tan,I,tol);

I = interval(1);
I_tan = tan(I);
I_true = interval(1.55740772465491);
res(end+1,1) = isequal(I_tan,I_true,tol);

% bounded, non-degenerate
I = interval(-1,0);
I_tan = tan(I);
I_true = interval(-1.55740772465491,0);
res(end+1,1) = isequal(I_tan,I_true,tol);

I = interval(-2,0);
I_tan = tan(I);
I_true = interval(-Inf,Inf);
res(end+1,1) = isequal(I_tan,I_true,tol);

I = interval(-pi/2,0);
I_tan = tan(I);
I_true = interval(-1.63312393531954e+16,0);
res(end+1,1) = isequal(I_tan,I_true,tol);

I = interval(-pi/2,pi/2);
I_tan = tan(I);
I_true = interval(-Inf,Inf);
res(end+1,1) = isequal(I_tan,I_true,tol);

% unbounded, non-degenerate
I = interval(-Inf,0);
I_tan = tan(I);
I_true = interval(-Inf,Inf);
res(end+1,1) = isequal(I_tan,I_true,tol);

I = interval(-Inf,Inf);
I_tan = tan(I);
I_true = interval(-Inf,Inf);
res(end+1,1) = isequal(I_tan,I_true,tol);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
