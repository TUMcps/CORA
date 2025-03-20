function res = test_emptySet_polytope
% test_emptySet_polytope - unit test function of polytope
%
% Syntax:
%    res = test_emptySet_polytope
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
% Written:       25-February-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

for n=1:10
    O = emptySet(n);
    P = polytope(O);
    assertLoop(representsa(P,'emptySet'),n)
end

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
