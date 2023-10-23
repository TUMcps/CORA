function res = isBigger(obj,E)
% isBigger - checks if ellipsoid(E2.Q) \subseteq ellipsoid(E1.Q) 
%
% Syntax:
%    res = isBigger(obj,E)
%
% Inputs:
%    obj - ellipsoid object
%    E - ellipsoid object
%
% Outputs:
%    res - boolean 
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
inputArgsCheck({{obj,'att','ellipsoid','scalar'};
                {E,'att','ellipsoid','scalar'}});

% check dimensions
if dim(obj) ~= dim(E)
    throw(CORAerror('CORA:dimensionMismatch',obj,E));
end

tol = min(obj.TOL,E.TOL);

% catch 1D case
if dim(obj)==1
    res = true;
    return;
end

% for simultaneous diagonalization, it is required that at least the first
% argument is pd (all eigenvalues >E.TOL)
% if that is not the case, then both are degenerate and either they share
% the same subspace so we can transform both into a subspace where at least
% the first is non-degenerate, or obj.Q cannot be bigger than E.Q

if ~isFullDim(obj) && isFullDim(E)
    res = false;
    return;
elseif ~isFullDim(obj) && ~isFullDim(E)
    % check if common subspace exists
    % 1/2 just for good measure
    Q_sum = 1/2*(obj.Q + E.Q);
    [U,S] = svd(Q_sum);

    r = rank(Q_sum);
    % if S does not contain at least 1 approx 0 eigenvalue, they do not
    % share common subspace
    if r==dim(obj)
        res = false;
        return;
    end
    
    % share common subspace
    % transform and cut
    obj = project(U'*obj,1:r);
    E = project(U'*E,1:r);
    
    
    % call isBigger again
    res = isBigger(obj,E);
    return;
end

% simulatenous diagonalization: Find Tb such that
% Tb'*Q1*Tb = I and Tb'*Q2*Tb = D (diagonal)
% if max(diag(D))<=1 => contained
[~,D] = simdiag(obj.Q,E.Q,tol);
tmp = max(diag(D));
res = tmp < 1+tol | withinTol(tmp,1+tol);

% ------------------------------ END OF CODE ------------------------------
