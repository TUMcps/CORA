function res = test_halfspace_uplus
% test_halfspace_uplus - unit test function of uplus
%
% Syntax:
%    res = test_halfspace_uplus
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

% plus
phs = +hs;
resvec(end+1) = all(phs.c == c, 'all');
resvec(end+1) = all(phs.d == d, 'all');

% compare with hs
resvec(end+1) = isequal(phs, hs);

% test empty case
resvec(end+1) = ~isemptyobject(+halfspace.empty(2));

% add results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
