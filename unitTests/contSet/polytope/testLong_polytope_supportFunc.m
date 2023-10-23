function res = testLong_polytope_supportFunc
% testLong_polytope_supportFunc - unit test function of supportFunc
%
% Syntax:
%    res = testLong_polytope_supportFunc
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

% tolerance
tol = 1e-8;

% number of tests
nrTests = 25;

% number of random directions
nrDir = 10;

for i=1:nrTests

    % random dimension
    n = randi(5);

    % instantiate random parallelotope
    Z = zonotope.generateRandom('Dimension',n,'NrGenerators',n);
    
    % convert to polytope
    P = polytope(Z);
    
    % compute support function for both along random directions and compare
    for j=1:nrDir
        % sample direction
        l = randn(n,1);

        % compute support function for zonotope
        sF_Z = supportFunc(Z,l);

        % compute support function for polytope
        sF_P = supportFunc(P,l);

        % compare
        if ~withinTol(sF_Z,sF_P,tol)
            res = false;
            return
        end
    end

    % instantiate random polytope
    P = polytope.generateRandom('Dimension',n);
    P = compact(P);
    
    % compute support function along all halfspaces
    nrCon = size(P.A,1);

    for j=1:nrCon
        % direction
        l = P.A(j,:)';

        % compute support function
        sF = supportFunc(P,l);

        % compare
        if ~withinTol(sF,P.b(j),tol)
            res = false;
            return
        end
    end


    % instantiate random polytope
    P = polytope.generateRandom('Dimension',n);
    P = compact(P);
    
    % add halfspace to make it degenerate
    randIdx = randi(size(P.A,1));
    dir = P.A(randIdx,:);
    val = P.b(randIdx);
    A = [P.A; -dir];
    b = [P.b; -val];
    P_ = polytope(A,b);

    if ~withinTol(supportFunc(P_,dir,'upper'),val) ...
            || ~withinTol(supportFunc(P_,dir,'lower'),val)
        res = false;
        return
    end
    

end

% ------------------------------ END OF CODE ------------------------------
