function res = test_nn_nnHelper_gerDerInterval()
% test_nn_nnHelper_gerDerInterval - tests the 
%     nnHelper.getDerInterval function
%
% Syntax:
%    res = test_nn_nnHelper_gerDerInterval()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnHelper/getDerInterval

% Authors:       Tobias Ladner
% Written:       17-February-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

l = -1;
u = 3;

% linear
coeffs = [2, 3];
[derl,deru] = nnHelper.getDerInterval(coeffs, l, u);
res = res && isequal(interval(derl,deru), interval(2, 2));

coeffs = [-5, 1];
[derl,deru] = nnHelper.getDerInterval(coeffs, l, u);
res = res && isequal(interval(derl,deru), interval(-5, -5));

% quadratic
coeffs = [1 0 2];
[derl,deru] = nnHelper.getDerInterval(coeffs, l, u);
res = res && isequal(interval(derl,deru), interval(-2, 6));

coeffs = [1 6 2];
[derl,deru] = nnHelper.getDerInterval(coeffs, l, u);
res = res && isequal(interval(derl,deru), interval(4, 12));

coeffs = [-1 7 2];
[derl,deru] = nnHelper.getDerInterval(coeffs, l, u);
res = res && isequal(interval(derl,deru), interval(1, 9));

% cubic
coeffs = [1 0.1 0 2];
[derl,deru] = nnHelper.getDerInterval(coeffs, l, u);
res = res && isequal(interval(derl,deru), interval(-1/300, 27.599999999999998));

end

% ------------------------------ END OF CODE ------------------------------
