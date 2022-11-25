function res_rand = testLongDuration_conZonotope_conZonotope
% testLongDuration_conZonotope_conZonotope - unit test function of
%    conZonotope (constructor)
%
% Syntax:  
%    res = testLongDuration_conZonotope_conZonotope
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

% Author:       Mark Wetzlinger
% Written:      19-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-12;
res_rand = true;
nrOfTests = 1000;
for i=1:nrOfTests

    % random dimension
    n = randi(25);
    % number of generators
    nrGens = randi(25);
    % number of constraints
    nrCon = randi(10);
    
    % random center, generator matrix (-> full matrix of zonotope),
    % as well as constraint matrix, constraint vector
    c = randn(n,1);
    G = randn(n,nrGens);
    Z = [c,G];
    A = randn(nrCon,nrGens);
    b = randn(nrCon,1);
    
    % admissible initializations
    % only zonotope matrix
    conZono = conZonotope(Z);
    if any(any(abs(conZono.Z - Z) > tol))
        res_rand = false; break;
    end

    % center and generator matrix
    conZono = conZonotope(c,G);
    if any(any(abs(conZono.Z - Z) > tol))
        res_rand = false; break;
    end

    % zonotope matrix, constraint matrix, constraint vector
    conZono = conZonotope(Z,A,b);
    if any(any(abs(conZono.Z - Z) > tol)) ...
            || any(any(abs(conZono.A - A) > tol)) ...
            || any(abs(conZono.b - b) > tol)
        res_rand = false; break;
    end
    
    % center, generator matrix, constraint matrix, constraint vector
    conZono = conZonotope(c,G,A,b);
    if any(any(abs(conZono.Z - Z) > tol)) ...
            || any(any(abs(conZono.A - A) > tol)) ...
            || any(abs(conZono.b - b) > tol)
        res_rand = false; break;
    end
    
    
    % wrong initializations
    c_plus1 = randn(n+1,1);
    c_mat = randn(n);
    G_plus1 = randn(n,nrGens+1);
    Z_plus1 = [c,G_plus1];
    A_plus1 = randn(nrCon,nrGens+1);
    b_plus1 = randn(nrCon+1,1);
    
    % center and generator matrix of different dimensions
    try
        conZono = conZonotope(c_plus1,G); % <- should throw error here
        res_rand = false; break;
    end
    
    % A does not fit Z
    try
        conZono = conZonotope(Z,A_plus1,b); % <- should throw error here
        res_rand = false; break;
    end
    try
        conZono = conZonotope(Z_plus1,A,b); % <- should throw error here
        res_rand = false; break;
    end
    
    % A does not fit b
    try
        conZono = conZonotope(Z,A,b_plus1); % <- should throw error here
        res_rand = false; break;
    end
    
    % center is a matrix
    if n ~= 1
        try
            conZono = conZonotope(c_mat,G); % <- should throw error here
            res_rand = false; break;
        end
    end
    
    % too many input arguments
    try
        conZono = conZonotope(c,G,A,b,b); % <- should throw error here
        res_rand = false; break;
    end
end


% combine results

if res_rand
    disp('testLongDuration_conZonotope_conZonotope successful');
else
    disp('testLongDuration_conZonotope_conZonotope failed');
end

%------------- END OF CODE --------------