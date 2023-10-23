function res = testLong_polytope_and
% testLong_polytope_and - unit test function of and
%
% Syntax:
%    res = testLong_polytope_and
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
% Written:       05-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% tolerance
tol = 1e-6;

% number of tests
nrTests = 25;

% 1D-4D intersection ------------------------------------------------------
% since we use the Minkowski sum, this test is restricted to low dimensions

for i=1:nrTests

    % random dimension
    n = randi(4);
    % random number of constraints
    nrCon = randi([2*n,2*n+2]);

    % instantiate polytope
    P1 = polytope.generateRandom('Dimension',n,'nrConstraints',nrCon);
    
    % compute shift that is required so that new polytope does not
    % intersect original one
    shift = 1.1*2*rad(interval(P1));
    P2 = P1 + shift;

    % compute intersection -> has to be empty
    P = and(P1,P2);

    % box has to contain the polytope
    if ~representsa(P, 'emptySet')
        res = false; return
    end

end


% nD intersection ---------------------------------------------------------

for i=1:nrTests

    % random dimension
    n = randi(10);
    % random number of constraints
    nrCon = randi([2*n,2*n+2]);

    % instantiate polytope
    P1 = polytope.generateRandom('Dimension',n,'nrConstraints',nrCon);
    
    % enlarge offset, so that new polytope encloses original one
    P2 = polytope(P1.A,P1.b+2);

    % compute intersection -> since P1 \subset P2, the intersection is P1
    P = and(P1,P2);

    % visualization
%     figure; hold on;
%     plot(P1);
%     plot(P2,[1,2],'r')
%     plot(P,[1,2],'g--');
%     close;

    % intersection has to be P1
    if ~eq(P1,P,tol)
        res = false; return
    end


    % starting from dimension 2
    n = randi([2,10]);
    % init row-normalized constraint matrices
    nrCon1 = randi([n+2,2*n+2]);
    A1 = randn(nrCon1,n);
    A1 = (A1' ./ vecnorm(A1'))';
    nrCon2 = randi([n+2,2*n+2]);
    A2 = randn(nrCon2,n);
    A2 = (A2' ./ vecnorm(A2'))';
    
    % same offset
    randVal = rand;
    b1 = randVal * ones(nrCon1,1);
    b2 = randVal * ones(nrCon2,1);

    % init polytopes
    P1 = polytope(A1,b1);
    P2 = polytope(A2,b2);

    % compute intersection
    P = P1 & P2;

    % true intersection
    P_ = polytope([A1; A2],[b1; b2]);

    % visualization
%     figure; hold on;
%     plot(P1);
%     plot(P2,[1,2],'k');
%     plot(P_,[1,2],'g');
%     plot(P,[1,2],'r--');
%     close;

    if ~eq(P,P_,tol)
        res = false; return
    end

end

% ------------------------------ END OF CODE ------------------------------
