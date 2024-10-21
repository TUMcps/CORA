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
    assertLoop(all(all(withinTol(E.Q,Q,tol))),i)

    % shape matrix and center
    E = ellipsoid(Q,q);
    assertLoop(all(all(withinTol(E.Q,Q,tol))),i)
    assertLoop(all(withinTol(E.q,q,tol)),i)
    
    
    Q_plus1 = temp * temp';
    
    % shape matrix non-psd (only n > 1)
    if n > 1
        assertThrowsAs(@ellipsoid,'CORA:wrongInputInConstructor',Q_nonpsd);
    end
    
    % shape matrix and center of different dimensions
    assertThrowsAs(@ellipsoid,'CORA:wrongInputInConstructor',Q,q_plus1);
    assertThrowsAs(@ellipsoid,'CORA:wrongInputInConstructor',Q_plus1,q);
    
    % center is a matrix
    if n > 1
        assertThrowsAs(@ellipsoid,'CORA:wrongValue',Q,q_mat);
    end
    
    % too many input arguments
    assertThrowsAs(@ellipsoid,'CORA:numInputArgsConstructor',Q,q,eps,q);
    
end

% ------------------------------ END OF CODE ------------------------------
