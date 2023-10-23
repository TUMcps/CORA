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

        % generate two zonotopes that intersect
        Z1 = zonotope.generateRandom('Dimension',n);
        Z2 = zonotope.generateRandom('Dimension',n);

        dir = rand(n,1);

        [~,x1] = supportFunc(Z1,dir);
        [~,x2] = supportFunc(Z2,-dir);

        tmp = center(Z1) - x1;
        x1_ = x1 + tmp/100;

        Z2 = Z2 + (x1_-x2);

        % check if GJK algorithm gives the correct result
        if ~GJKalgorithm(Z1,Z2)
            throw(CORAerror('CORA:testFailed'));
        end

        % generate two zonotopes that do not intersect
        Z1 = zonotope.generateRandom('Dimension',n);
        Z2 = zonotope.generateRandom('Dimension',n);

        dir = rand(n,1);

        [~,x1] = supportFunc(Z1,dir);
        [~,x2] = supportFunc(Z2,-dir);

        tmp = x1 - center(Z1);
        x1_ = x1 + tmp/100;

        Z2 = Z2 + (x1_-x2);

        % check if GJK algorithm gives the correct result
        if GJKalgorithm(Z1,Z2)
            throw(CORAerror('CORA:testFailed'));
        end
    end
end

res = true;

% ------------------------------ END OF CODE ------------------------------
