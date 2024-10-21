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
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann, Mark Wetzlinger
% Written:       26-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% empty ellipsoid
E = ellipsoid.empty(2);
assert(representsa_(E,'emptySet',eps))

tol = 1e-12;

% init cases
E1 = ellipsoid([ 5.4387811500952807 12.4977183618314545 ; 12.4977183618314545 29.6662117284481646 ], [ -0.7445068341257537 ; 3.5800647524843665 ], 0.000001);

Q = E1.Q;
q = center(E1);
n = length(q);

% only shape matrix
E = ellipsoid(Q);
assert(all(all(withinTol(E.Q,Q,tol))))

% shape matrix and center
E = ellipsoid(Q,q);
assert(all(all(withinTol(E.Q,Q,tol))))
assert(all(all(withinTol(E.q,q,tol))))    

% wrong instantiations
% shape matrix non-psd (only n > 1)
if n > 1
    % random non-psd matrix
    Q_nonpsd = randn(n);
    [U,S,V] = svd(Q_nonpsd);
    ind_p = diag(S)>0;
    if sum(ind_p)==n
        i_r = randi([1,n]);
        S(i_r,i_r) = -S(i_r,i_r);
    end
    assertThrowsAs(@ellipsoid,'CORA:wrongInputInConstructor',U*S*V');
end

% shape matrix and center of different dimensions
assertThrowsAs(@ellipsoid,'CORA:wrongInputInConstructor',Q,[q;1]);
assertThrowsAs(@ellipsoid,'CORA:wrongInputInConstructor',blkdiag(Q,1),q);

% center is a matrix
if n ~= 1
    assertThrowsAs(@ellipsoid,'CORA:wrongValue',Q,repmat(q,1,2));
end

% too many input arguments
assertThrowsAs(@ellipsoid,'CORA:numInputArgsConstructor',Q,q,eps,q);
    
end

% ------------------------------ END OF CODE ------------------------------
