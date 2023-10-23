function res = testLong_polytope_plus
% testLong_polytope_plus - unit test function of plus
%
% Syntax:
%    res = testLong_polytope_plus
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

% Authors:       Mark Wetzlinger
% Written:       01-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 25;

for i=1:nrTests

    % random dimension
    n = randi(4);

    % instantiate polytope
    P = polytope.generateRandom('Dimension',n);
    
    % choose random constraint vector
    randCon = randi(size(P.A,1));
    con = P.A(randCon,:)';

    % normalize to length 1
    origLength = vecnorm(con);
    con = con ./ origLength;

    % shift polytope by this vector
    P_ = P + con;

    % compare offsets
    if ~withinTol(P.b(randCon)+origLength,P_.b(randCon),1e-10)
        res = false; break
    end

end

% ------------------------------ END OF CODE ------------------------------
