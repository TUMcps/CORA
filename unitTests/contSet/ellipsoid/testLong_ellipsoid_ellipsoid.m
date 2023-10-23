function res = testLong_ellipsoid_ellipsoid
% testLong_ellipsoid_ellipsoid - unit test function of ellipsoid
%
% Syntax:
%    res = testLong_ellipsoid_ellipsoid
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
nrOfTests = 100;
for i=1:nrOfTests
    %%% generate all variables necessary to replicate results
    % random dimension
    n = randi(15);
    % random shape matrix (psd and non-psd) and random center
    q = randn(n,1);
    Q_nonpsd = randn(n);
    % wrong initializations
    q_plus1 = randn(n+1,1);
    q_mat = randn(n);
    temp = randn(n+1);
    %%%
    Q = Q_nonpsd * Q_nonpsd';
    
    % admissible initializations
    % only shape matrix
    E = ellipsoid(Q);
    if ~all(all(withinTol(E.Q,Q,tol)))
        res = false; break;
    end

    % shape matrix and center
    E = ellipsoid(Q,q);
    if ~all(all(withinTol(E.Q,Q,tol))) || ~all(withinTol(E.q,q,tol))
        res = false; break;
    end
    
    
    Q_plus1 = temp * temp';
    
    % shape matrix non-psd (only n > 1)
    if n > 1
        try
            E = ellipsoid(Q_nonpsd); % <- should throw error here
            res = false; break;
        end
    end
    
    % shape matrix and center of different dimensions
    try
        E = ellipsoid(Q,q_plus1); % <- should throw error here
        res = false; break;
    end
    try
        E = ellipsoid(Q_plus1,q); % <- should throw error here
        res = false; break;
    end
    
    % center is a matrix
    if n ~= 1
        try
            E = ellipsoid(Q,q_mat); % <- should throw error here
            res = false; break;
        end
    end
    
    % too many input arguments
    try
        E = ellipsoid(Q,q,eps,q); % <- should throw error here
        res = false; break;
    end
    
end

% ------------------------------ END OF CODE ------------------------------
