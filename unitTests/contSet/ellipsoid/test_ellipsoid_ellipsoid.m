function res = test_ellipsoid_ellipsoid
% test_ellipsoid_ellipsoid - unit test function of ellipsoid
%
% Syntax:  
%    res = test_ellipsoid_ellipsoid
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

% Author:       Victor Gassmann, Mark Wetzlinger
% Written:      26-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
tol = 1e-12;
load cases.mat E_c
for i=1:length(E_c)
    E1 = E_c{i}.E1; % non-deg
    Q = E1.Q;
    q = center(E1);
    n = length(q);
    
    % only shape matrix
    E = ellipsoid(Q);
    if any(any(abs(E.Q - Q) > tol))
        res = false; break;
    end

    % shape matrix and center
    E = ellipsoid(Q,q);
    if any(any(abs(E.Q - Q) > tol)) || any(abs(E.q - q) > tol)
        res = false; break;
    end
    
    
    % shape matrix non-psd (only n > 1)
    if n > 1
        try
            % random non-psd matrix
            Q_nonpsd = randn(n);
            [U,S,V] = svd(Q_nonpsd);
            ind_p = diag(S)>0;
            if sum(ind_p)==n
                i_r = randi([1,n]);
                S(i_r,i_r) = -S(i_r,i_r);
            end
            E = ellipsoid(U*S*V'); % <- should throw error here
            res = false; break;
        end
    end
    
    % shape matrix and center of different dimensions
    try
        E = ellipsoid(Q,[q;1]); % <- should throw error here
        res = false; break;
    end
    try
        E = ellipsoid(blkdiag(Q,1),q); % <- should throw error here
        res = false; break;
    end
    
    % center is a matrix
    if n ~= 1
        try
            E = ellipsoid(Q,repmat(q,1,2)); % <- should throw error here
            res = false; break;
        end
    end
    
    % too many input arguments
    try
        E = ellipsoid(Q,q,eps,q); % <- should throw error here
        res = false; break;
    end
    
    
end


if res
    disp([mfilename,' successful']);
else
    disp([mfilename,' failed']);
end
%------------- END OF CODE --------------
