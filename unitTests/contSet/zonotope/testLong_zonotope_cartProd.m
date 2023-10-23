function res = testLong_zonotope_cartProd
% testLong_zonotope_cartProd - unit test function of Cartesian
%    product
%
% Syntax:
%    res = testLong_zonotope_cartProd
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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       03-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% number of tests
nrTests = 1000;
res = true;

for i=1:nrTests
    
    % zonotope-zonotope case:

    % random dimensions
    n1 = randi(10);
    n2 = randi(10);

    % random number of generators
    m1 = randi(n1*2);
    m2 = randi(n2*2);

    % random centers and generator matrices
    c1 = randn(n1,1);
    c2 = randn(n2,1);
    G1 = randn(n1,m1);
    G2 = randn(n2,m2);

    % instantiate zonotopes
    Z1 = zonotope(c1,G1);
    Z2 = zonotope(c2,G2);
    
    % compute Cartesian product
    Z_ = cartProd(Z1,Z2);
    
    % obtain center and generator matrix
    c = center(Z_);
    G = generators(Z_);
    
    % check result
    if ~( compareMatrices(c,[c1;c2]) && compareMatrices(G,blkdiag(G1,G2)) )
        res = false; break
    end


    % zonotope-interval case

    % random dimensions
    n1 = randi(10);
    n2 = randi(10);

    % random number of generators
    m1 = randi(n1*2);

    % random center and generator matrix
    c1 = randn(n1,1);
    G1 = randn(n1,m1);
    
    % random lower and upper bounds
    lb2 = -rand(n2,1);
    ub2 = rand(n2,1);

    % instantiate zonotope and interval
    Z1 = zonotope(c1,G1);
    I2 = interval(lb2,ub2);
    
    % compute Cartesian product
    Z_ = cartProd(Z1,I2);
    
    % obtain center and generator matrix
    c = center(Z_);
    G = generators(Z_);
    
    % check result
    if ~( compareMatrices(c,[c1;0.5*(lb2+ub2)]) ...
            && compareMatrices(G,blkdiag(G1,0.5*diag(ub2-lb2))) )
        res = false; break
    end

    % zonotope-numeric case

    % random dimensions
    n1 = randi(10);
    n2 = randi(10);

    % random number of generators
    m1 = randi(n1*2);

    % random center and generator matrix
    c1 = randn(n1,1);
    G1 = randn(n1,m1);
    Z1 = zonotope(c1,G1);
    
    % random numeric vector
    num = randn(n2,1);
    
    % compute Cartesian product
    Z_ = cartProd(Z1,num);
    
    % obtain center and generator matrix
    c = center(Z_);
    G = generators(Z_);
    
    % check result
    if ~( compareMatrices(c,[c1;num]) ...
            && compareMatrices(G,[G1;zeros(n2,m1)]) )
        res = false; break
    end

    % compute Cartesian product
    Z_ = cartProd(num,Z1);
    
    % obtain center and generator matrix
    c = center(Z_);
    G = generators(Z_);
    
    % check result
    if ~( compareMatrices(c,[num;c1]) ...
            && compareMatrices(G,[zeros(n2,m1);G1]) )
        res = false; break
    end

end

% ------------------------------ END OF CODE ------------------------------
