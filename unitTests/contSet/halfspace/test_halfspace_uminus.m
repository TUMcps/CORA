function res = test_halfspace_uminus
% test_halfspace_uminus - unit test function of uminus
%
% Syntax:
%    res = test_halfspace_uminus
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
% See also: -

% Authors:       Tobias Ladner
% Written:       06-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = true(0);

% init
c = [1 1];
d = 2;
hs = halfspace(c,d);

% negate
nhs = -hs;
resvec(end+1) = all(nhs.c == -c, 'all');
resvec(end+1) = all(nhs.d == d, 'all');

% compare with -1 * hs
resvec(end+1) = isequal(nhs, -1*hs);

% test empty case
resvec(end+1) = ~isemptyobject(-halfspace.empty(2));

% add results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
