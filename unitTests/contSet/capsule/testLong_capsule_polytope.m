function res = testLong_capsule_polytope
% testLong_capsule_polytope - unit test function of polytope conversion
%
% Syntax:
%    res = testLong_capsule_polytope
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
% Written:       25-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% tolerance
tol = 1e-10;

% number of tests
nrTests = 50;

for j=1:nrTests

    % random dimension
    n = randi(10);
    
    % random capsule
    C = capsule.generateRandom("Dimension",n);

    % convert to polytope
    P = polytope(C,'outer');

    % check support function evaluation in 2n axis-aligned directions
    for i=1:n
        for e=[1,-1]
            dir = [zeros(i-1,1); e; zeros(n-i,1)];
            val_C = supportFunc_(C,dir,'upper');
            val_P = supportFunc_(P,dir,'upper');
            if val_P < val_C
                assertLoop(withinTol(val_P,val_C,tol),j,i,e)
            end
        end
    end

    % capsule with only center
    C_center = capsule(C.c);
    P_center = polytope(C_center);
    assert(contains(P_center,C_center.c,'exact',tol))
%   assert(representsa_(P_center,'point',tol))

    C_generator = capsule(C.c,C.g);
    P_generator = polytope(C_generator);
    % check random support function evaluation
    dir = randn(n,1);
    assert(withinTol(supportFunc_(C_generator,dir,'upper'),...
            supportFunc_(P_generator,dir,'upper'),tol))

end

end

% ------------------------------ END OF CODE ------------------------------
