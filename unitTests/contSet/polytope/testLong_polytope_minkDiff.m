function res = testLong_polytope_minkDiff
% testLong_polytope_minkDiff - unit test function of minkDiff
%
% Syntax:
%    res = testLong_polytope_minkDiff
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
    n = randi(5);

    % instantiate polytope
    P2 = polytope.generateRandom('Dimension',n);
    
    % shift offset vector by a little
    z = min(abs(P2.b));
    b_ = P2.b + z;

    % instantiate larger polytope, enclosing P2
    P1 = polytope(P2.A,b_);

    % compute Minkowski difference
    P = minkDiff(P1,P2);

    % visualization
%     figure; hold on;
%     h1 = plot(P1);
%     h2 = plot(P2,[1,2],'k');
%     hdiff = plot(P,[1,2],'r');
%     legend([h1,h2,hdiff],'P1','P2','minkDiff');

    % result has to contain the origin as P2 \subset P1
    if ~contains(P,zeros(n,1))
        throw(CORAerror('CORA:testFailed'));
    end

    % compute Minkowski difference the other way around
    P = minkDiff(P2,P1);

    % result has to be empty as P2 \subset P1
    if ~representsa(P, 'emptySet')
        throw(CORAerror('CORA:testFailed'));
    end


    % subtract an interval or a zonotope from a polytope
    P1 = polytope.generateRandom('Dimension',n,'NrConstraints',2*n);
    I = interval.generateRandom('Dimension',n);
    Z = zonotope.generateRandom('Dimension',n,'Center',zeros(n,1));

    % should just run through...
    P = minkDiff(P1,I);
    P = minkDiff(P1,Z);

end

% ------------------------------ END OF CODE ------------------------------
