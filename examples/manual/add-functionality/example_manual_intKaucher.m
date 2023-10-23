function example_manual_intKaucher()
% example_manual_intKaucher - example from the manual demontrating the
% intKaucher operation as defined in the manual
%
% Syntax:
%   example_manual_intKaucher()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Tobias Ladner
% Written:       27-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% function f
f = @(x) x^2 - x;

% compute gradient
syms x;
df = gradient(f(x));
df = matlabFunction(df);

% compute bounds for gradient
I = interval(2,3);
c = center(I);
gr = df(I);

% compute inner-approximation of the range
x = intKaucher(3,2);
gr = intKaucher(infimum(gr), supremum(gr));

res = f(c) + gr*(x - c);

% ------------------------------ END OF CODE ------------------------------
