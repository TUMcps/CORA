function res = isBigger(E1,E2)
% isBigger - checks if an ellipsoid is bigger than another ellipsoid when
%    both centers are moved to the origin, i.e.,
%    ellipsoid(E2.Q) \subseteq ellipsoid(E1.Q) 
%
% Syntax:
%    res = isBigger(E1,E2)
%
% Inputs:
%    E1 - ellipsoid object
%    E2 - ellipsoid object
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       10-June-2022
% Last update:   20-March-2023 (VG, allow degeneracy)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input
inputArgsCheck({{E1,'att','ellipsoid','scalar'};
                {E2,'att','ellipsoid','scalar'}});

% check dimensions
equalDimCheck(E1,E2);

% set tolerance
tol = min(E1.TOL,E2.TOL);

% catch 1D case
if dim(E1)==1
    res = true;
    return;
end

% for simultaneous diagonalization, it is required that at least the first
% argument is positive definite (i.e., all eigenvalues >E.TOL)
% if that is not the case, then both are degenerate and either they share
% the same subspace so we can transform both into a subspace where at least
% the first is non-degenerate, or E1.Q cannot be bigger than E2.Q

if ~isFullDim(E1) && isFullDim(E2)
    res = false;
    return;

elseif ~isFullDim(E1) && ~isFullDim(E2)
    % check if common subspace exists (factor 1/2 just for good measure)
    Q_sum = 1/2*(E1.Q + E2.Q);
    [U,S,~] = svd(Q_sum);

    r = rank(Q_sum);
    % if S does not contain at least 1 approx. 0 eigenvalue, they do not
    % share common subspace
    if r==dim(E1)
        res = false;
        return;
    end
    
    % share common subspace
    % transform and cut
    E1 = project(U'*E1,1:r);
    E2 = project(U'*E2,1:r);
    
    % call isBigger again
    res = isBigger(E1,E2);
    return;
end

% simultaneous diagonalization: Find Tb such that
%   Tb'*Q1*Tb = I and Tb'*Q2*Tb = D (diagonal)
[~,D] = simdiag(E1.Q,E2.Q,tol);
tmp = max(diag(D));

% if max(diag(D)) <= 1 => contained
res = tmp < 1+tol | withinTol(tmp,1+tol);

% ------------------------------ END OF CODE ------------------------------
