function res = testLong_capsule_capsule
% testLong_capsule_capsule - unit test function of capsule (constructor)
%
% Syntax:
%    res = testLong_capsule_capsule
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

tol = 1e-12;

res = true;
nrOfTests = 1000;
for i=1:nrOfTests

    % random dimension
    n = randi([2,25]);
    
    % random center, generator, and radius
    c = randn(n,1);
    g = randn(n,1);
    r = rand(1);
    
    % admissible initializations
    % only center
    C = capsule(c);
    if ~compareMatrices(C.c,c)
        res = false; break;
    end

    % center and generator
    C = capsule(c,g);
    if ~compareMatrices(C.c,c) || ~compareMatrices(C.g,g)
        res = false; break;
    end
    
    % center, generator, and radius
    C = capsule(c,g,r);
    if ~compareMatrices(C.c,c) || ~compareMatrices(C.g,g) || ~withinTol(C.r,r)
        res = false; break;
    end
    
    % prevent 1D cases (where g can be interpreted as r)
    n = n + 1;
    c = randn(n,1);
    g = randn(n,1);
    
    % wrong initializations
    cn_1 = randn(n+1,1);
    gn_1 = randn(n+1,1);
    rneg = -rand(1);
    rvec = rand(n,1);
    
    % mismatch between center and generator
    try
        C = capsule(c,gn_1); % <- should throw error here
        res = false; break;
    end
    try
        C = capsule(cn_1,g); % <- should throw error here
        res = false; break;
    end
    
    % negative radius
    try
        C = capsule(c,g,rneg); % <- should throw error here
        res = false; break;
    end
    
    % radius as a vector
    try
        C = capsule(c,g,rvec); % <- should throw error here
        res = false; break;
    end
    
    % too many input arguments
    try
        C = capsule(c,g,r,r); % <- should throw error here
        res = false; break;
    end
    
end

% ------------------------------ END OF CODE ------------------------------
