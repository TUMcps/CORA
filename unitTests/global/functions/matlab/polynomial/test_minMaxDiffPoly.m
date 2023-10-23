function res = test_minMaxDiffPoly()
% test_minMaxDiffPoly - tests the minMaxDiffPoly function
%
% Syntax:
%    res = test_minMaxDiffPoly()
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
% See also: -

% Authors:       Tobias Ladner
% Written:       30-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% simple cases

[diffl,diffu] = minMaxDiffPoly([1 0],[0], -1, 2);
L = interval(diffl,diffu);
resvec(end+1) = isequal(L,interval(-1,2));

[diffl,diffu] = minMaxDiffPoly([0 0],[1 0], -1, 2);
L = interval(diffl,diffu);
resvec(end+1) = isequal(L,interval(-2,1));

[diffl,diffu] = minMaxDiffPoly([2 0 0],[1 0], -1, 1);
L = interval(diffl,diffu);
resvec(end+1) = isequal(L,interval(-0.125,3));

[diffl,diffu] = minMaxDiffPoly([2 3],[2 0], -1, 1);
L = interval(diffl,diffu);
resvec(end+1) = isequal(L,interval(3,3));

% higher-order polynomials

[diffl,diffu] = minMaxDiffPoly([4 2 1 6],[-2 4], -3, 1);
L = interval(diffl,diffu);
resvec(end+1) = isequal(L,interval(-97,11));

[diffl,diffu] = minMaxDiffPoly([-1 2 -1 6],[2 5 4], 2, 5);
L = interval(diffl,diffu);
resvec(end+1) = isequal(L,interval(-153,-18));

% test swapped
coeffs1 = [2 -9 3];
coeffs2 = [-5 4 6 2];
l = -6;
u = -4;
[diffl1,diffu1] = minMaxDiffPoly(coeffs1,coeffs2,l,u);
L1 = interval(diffl1,diffu1);
[diffl2,diffu2] = minMaxDiffPoly(coeffs2,coeffs1,l,u);
L2 = interval(diffl2,diffu2);
resvec(end+1) = isequal(L1,-L2);

% gather results
res = all(resvec);

end

% ------------------------------ END OF CODE ------------------------------
