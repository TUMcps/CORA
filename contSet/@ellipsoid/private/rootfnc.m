function [y,Q,q] = rootfnc(p,W1,q1,W2,q2)
% rootfnc - polynomial whose root corresponds to the correct p (see 'and_')
%
% Syntax:
%    [y,Q,q] = rootfnc(p,W1,q1,W2,q2)
%
% Inputs:
%    p - current p value
%    W1,q1,W2,q2 - shape matrices and centers of E1,E2 (see 'and')
%
% Outputs:
%    y - current function value (should go to 0)
%    Q,q - current solution, i.e., ellipsoid(Q,q)
%
% References: Very heavily inspired by file 'ell_fusionlambda.m' of the
%    Ellipsoidal toolbox:
%    https://de.mathworks.com/matlabcentral/fileexchange/21936-ellipsoidal-toolbox-et
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ellipsoid/and_

% Authors:       Victor Gassmann
% Written:       14-October-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

n = size(W1,1);
X = p*W1+(1-p)*W2;
X_inv = inv(X);
X_inv = 0.5*(X_inv+X_inv');
a = 1-p*(1-p)*(q2-q1)'*W2*X_inv*W1*(q2-q1);
q = X_inv*(p*W1*q1+(1-p)*W2*q2);
Q = a*inv(X);
detX = det(X);
%y = a*detX*trace(detX*X_inv*(W1 - W2)) - n*detX^2* ...
%      (2*q'*W1*q1 - 2*q'*W2*q2 + q'*(W2 - W1)*q - q1'*W1*q1 + q2'*W2*q2);
y = a*detX^2*trace(X_inv*(W1-W2)) - n*detX^2*(2*q'*(W1*q1-W2*q2) ...
    + q'*(W2-W1)*q - q1'*W1*q1+q2'*W2*q2);

% ------------------------------ END OF CODE ------------------------------
