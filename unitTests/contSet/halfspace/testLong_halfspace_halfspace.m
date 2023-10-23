function res = testLong_halfspace_halfspace
% testLong_halfspace_halfspace - unit test function of halfspace
%
% Syntax:
%    res = testLong_halfspace_halfspace
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
% Written:       20-March-2021
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
    
    % random normal vector, offset
    c = randn(n,1);
    d = randn(1);
    
    % admissible initialization
    hs = halfspace(c,d);
    if ~compareMatrices(hs.c,c) || ~withinTol(hs.d,d)
        res = false; break;
    end
    hs = halfspace(c',d);
    if ~compareMatrices(hs.c,c) || ~withinTol(hs.d,d)
        res = false; break;
    end
    
    % wrong initializations
    c_mat = randn(n+1);
    d_vec = randn(n+1,1);
    d_mat = randn(n+1);
    
    % only one input
    try
        hs = halfspace(c); % <- should throw error here
        res = false; break;
    end
    
    % normal vector is a matrix
    try
        hs = halfspace(c_mat,d); % <- should throw error here
        res = false; break;
    end
    
    % offset is a vector
    try
        hs = halfspace(c,d_vec); % <- should throw error here
        res = false; break;
    end
    
    % offset is a matrix
    try
        hs = halfspace(c,d_mat); % <- should throw error here
        res = false; break;
    end
    
    % too many input arguments
    try
        hs = halfspace(c,d,d); % <- should throw error here
        res = false; break;
    end 
end

% ------------------------------ END OF CODE ------------------------------
