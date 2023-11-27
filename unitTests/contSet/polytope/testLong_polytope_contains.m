function res = testLong_polytope_contains
% testLong_polytope_contains - unit test function for containment
%    check
%
% Syntax:
%    res = testLong_polytope_contains()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Viktor Kotsev
% Written:       ---
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
   
    % create random polytope
    P = polytope.generateRandom('Dimension',n);
    
    % check if randomly generated points are inside
    Y = randPoint(P,n);
    if ~contains(P,Y,'exact',1e-12)
        throw(CORAerror('CORA:testFailed'));
    end


    % instantiate polytope
    P2 = polytope.generateRandom('Dimension',n);
    b_ = P2.b + 1;
    P1 = polytope(P2.A,b_);
    
    % check containment
    if ~contains(P1,P2)
        throw(CORAerror('CORA:testFailed'));
    end
    
    % generate zonotope
    Z = zonotope.generateRandom('Dimension',n);
    E = ellipsoid(Z,'inner:norm');
    P = polytope(Z);
    % check if Ei is in P
    if ~contains(P,E,'exact',E.TOL)
        throw(CORAerror('CORA:testFailed'));
    end
end

% ------------------------------ END OF CODE ------------------------------
