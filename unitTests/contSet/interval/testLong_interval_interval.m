function res = testLong_interval_interval
% testLong_interval_interval - unit test function of interval
%
% Syntax:
%    res = testLong_interval_interval
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

tol = 1e-12;

res = true;
nrOfTests = 1000;
for i=1:nrOfTests

    % random dimension
    n = randi(25);
    
    % random lower bound, random upper bound
    a = -rand(n,1);
    b = rand(n,1);
    a_mat = -rand(n+1);
    b_mat = rand(n+1);
    
    % admissible initializations
    I = interval(a,b);
    if any(abs(I.inf - a) > tol) || any(abs(I.sup - b) > tol)
        res = false; break;
    end
    
    I = interval(a);
    if any(abs(I.inf - a) > tol) || any(abs(I.sup - a) > tol)
        res = false; break;
    end
    
    I = interval(a_mat);
    if any(any(abs(I.inf - a_mat) > tol)) || any(any(abs(I.sup - a_mat) > tol))
        res = false; break;
    end
    
    I = interval(a_mat,b_mat);
    if any(any(abs(I.inf - a_mat) > tol)) || any(any(abs(I.sup - b_mat) > tol))
        res = false; break;
    end
    
    % wrong initializations
    a_large = 1+rand(n,1);
    b_small = -1-rand(n,1);
    a_plus1 = -rand(n+1,1);
    b_plus1 = rand(n+1,1);
    
    % lower limit larger than upper limit
    try
        I = interval(a,b_small); % <- should throw error here
        res = false; break;
    end
    try
        I = interval(a_large,b); % <- should throw error here
        res = false; break;
    end
    
    % size of limits do not match
    try
        I = interval(a_plus1,b); % <- should throw error here
        res = false; break;
    end
    try
        I = interval(a,b_plus1); % <- should throw error here
        res = false; break;
    end
    try
        I = interval(a_mat,b); % <- should throw error here
        res = false; break;
    end
    try
        I = interval(a,b_mat); % <- should throw error here
        res = false; break;
    end
    
    % too many input arguments
    try
        I = interval(a,b,b); % <- should throw error here
        res = false; break;
    end 
end

% ------------------------------ END OF CODE ------------------------------
