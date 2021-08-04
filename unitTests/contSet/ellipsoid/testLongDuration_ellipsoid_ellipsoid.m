function res = testLongDuration_ellipsoid_ellipsoid
% testLongDuration_ellipsoid_ellipsoid - unit test function of ellipsoid
%
% Syntax:  
%    res = testLongDuration_ellipsoid_ellipsoid
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

% empty conHyperplane
E = ellipsoid();
res_empty = true;
if ~isempty(E)
    res_empty = false;
end

res_rand = true;
nrOfTests = 100;
for i=1:nrOfTests

    % random dimension
    n = randi(15);
    
    % random shape matrix (psd and non-psd) and random center
    q = randn(n,1);
    Q_nonpsd = randn(n);
    Q = Q_nonpsd * Q_nonpsd';
    
    % admissible initializations
    % only shape matrix
    E = ellipsoid(Q);
    if any(any(abs(E.Q - Q) > tol))
        res_rand = false; break;
    end

    % shape matrix and center
    E = ellipsoid(Q,q);
    if any(any(abs(E.Q - Q) > tol)) || any(abs(E.q - q) > tol)
        res_rand = false; break;
    end
    
    
    % wrong initializations
    q_plus1 = randn(n+1,1);
    q_mat = randn(n);
    temp = randn(n+1);
    Q_plus1 = temp * temp';
    
    % shape matrix non-psd (only n > 1)
    if n > 1
        try
            E = ellipsoid(Q_nonpsd); % <- should throw error here
            res_rand = false; break;
        end
    end
    
    % shape matrix and center of different dimensions
    try
        E = ellipsoid(Q,q_plus1); % <- should throw error here
        res_rand = false; break;
    end
    try
        E = ellipsoid(Q_plus1,q); % <- should throw error here
        res_rand = false; break;
    end
    
    % center is a matrix
    if n ~= 1
        try
            E = ellipsoid(Q,q_mat); % <- should throw error here
            res_rand = false; break;
        end
    end
    
    % too many input arguments
    try
        E = ellipsoid(Q,q,eps,q); % <- should throw error here
        res_rand = false; break;
    end
    
end


% combine results
res = res_empty && res_rand;

if res
    disp('testLongDuration_ellipsoid_ellipsoid successful');
else
    disp('testLongDuration_ellipsoid_ellipsoid failed');
end

%------------- END OF CODE --------------