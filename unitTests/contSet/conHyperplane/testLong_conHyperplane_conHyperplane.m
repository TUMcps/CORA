function res = testLong_conHyperplane_conHyperplane
% testLong_conHyperplane_conHyperplane - unit test function of
%    conHyperplane (constructor)
%
% Syntax:
%    res = testLong_conHyperplane_conHyperplane
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
% Written:       19-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrOfTests = 1000;

for i=1:nrOfTests

    % random dimension
    n = randi(25);
    % number of constraints
    nrCon = randi(25);
    
    % random normal vector, offset, constraint matrix, constraint vector
    a = randn(1,n);
    b = randn(1);
    C = randn(nrCon,n);
    d = randn(nrCon,1);
    
    % admissible initializations

    % normal vector and offset
    hyp = conHyperplane(a,b);
    if ~all(withinTol(hyp.a,a)) || ~withinTol(hyp.b,b)
        throw(CORAerror('CORA:testFailed'));
    end

    % equality constraint and inequality constraints
    hyp = conHyperplane(a,b,C,d);
    if ~all(withinTol(hyp.a,a)) || ~withinTol(hyp.b,b) ...
            || ~compareMatrices(hyp.C,C) || ~compareMatrices(hyp.d,d)
        throw(CORAerror('CORA:testFailed'));
    end
        
    % normal vector, offset, and constraint matrix, constraint vector
    hyp = conHyperplane(a,b,C,d);
    if ~all(withinTol(hyp.a,a)) || ~withinTol(hyp.b,b) ...
            || ~compareMatrices(hyp.C,C) || ~compareMatrices(hyp.d,d)
        throw(CORAerror('CORA:testFailed'));
    end
    
    
    % wrong initializations
    a_plus1 = randn(n+1,1);
    b_vec = randn(n+1,1);
    C_plus1 = randn(nrCon,n+1);
    d_plus1 = randn(nrCon+1,1);
    
    % offset as vector
    try
        hyp = conHyperplane(a,b_vec); % <- should throw error here
        throw(CORAerror('CORA:testFailed'));
    end
    
    % a does not fit C
    try
        hyp = conHyperplane(a,b,C_plus1,d); % <- should throw error here
        throw(CORAerror('CORA:testFailed'));
    end 
    
    % a does not fit d
    try
        hyp = conHyperplane(a,b,C,d_plus1); % <- should throw error here
        throw(CORAerror('CORA:testFailed'));
    end 
    
    % too many input arguments
    try
        hyp = conHyperplane(a,b,C,d,d); % <- should throw error here
        throw(CORAerror('CORA:testFailed'));
    end 
end

% ------------------------------ END OF CODE ------------------------------
