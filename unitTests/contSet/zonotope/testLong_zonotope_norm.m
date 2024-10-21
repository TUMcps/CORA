function res = testLong_zonotope_norm
% testLong_zonotope_norm - unit test function of norm
%
% Syntax:
%    res = testLong_zonotope_norm
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

% Authors:       Victor Gassmann
% Written:       31-July-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

TOL = 1e-6;

% loop over dimensions
for i=2:4

    % loop over number of generators
    for j=i:2:10

        % instantiate random zonotope
        Z = zonotope([zeros(i,1),10*randn(i,j)]);

        % 2 norm test
        val2_exact = norm(Z,2,'exact');
        val2_ub = norm(Z,2,'ub');
        val2_ubc = norm(Z,2,'ub_convex');

        % compute vertices
        V = vertices(Z);

        % check exact vs. upper bound
        if val2_exact > val2_ub 
            assertLoop(withinTol(val2_exact,val2_ub,TOL),i,j)
        end
        
        % check exact vs. upper bound (convex)
        if val2_exact > val2_ubc 
            assertLoop(withinTol(val2_exact,val2_ubc,TOL),i,j)
        end
        
        % check exact vs. norm of all vertices
        val = abs(val2_exact-max(sqrt(sum(V.^2))))/val2_exact;
        if val > 0 
            assertLoop(withinTol(val,0,TOL),i,j)
        end

    end

end

% ------------------------------ END OF CODE ------------------------------
