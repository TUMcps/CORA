function res = testLong_GJKalgorithm()
% testLong_GJKalgorithm - unit test function for the GJK algorithm
%
% Syntax:
%    res = testLong_GJKalgorithm()
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
% See also: GJKalgorithm

% Authors:       Niklas Kochdumper
% Written:       21-April-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

dims = [2 3 4 10];

for n = dims    
    for i = 1:10
        disp("n: " + n + ", i = " + i);

        % generate two zonotopes: we compute the support vectors in
        % opposite directions and shift a zonotope by the vector between
        % them; this would make the zonotopes meet at exactly one point, so
        % we shift one support vector by a tiny bit to make them either
        % intersect or not intersect

        Z1 = zonotope.generateRandom('Dimension',n);
        Z2 = zonotope.generateRandom('Dimension',n);

        dir = randn(n,1);
        [~,x1] = supportFunc(Z1,dir);
        [~,x2] = supportFunc(Z2,-dir);
        offset = center(Z1) - x1;
        Z2_intersects = Z2 + (x1 + offset/100 - x2);
        Z2_notintersects = Z2 + (x1 - offset/100 - x2);

        if ~GJKalgorithm(Z1,Z2_intersects)
            throw(CORAerror('CORA:testFailed'));
        end

        if GJKalgorithm(Z1,Z2_notintersects)
            throw(CORAerror('CORA:testFailed'));
        end
    end
end

res = true;

% ------------------------------ END OF CODE ------------------------------
