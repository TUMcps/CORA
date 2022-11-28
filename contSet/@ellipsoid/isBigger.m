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

% Author:       Victor Gassmann
% Written:      10-June-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% check input
inputArgsCheck({{obj,'att','ellipsoid','scalar'};
                {E,'att','ellipsoid','scalar'}});

% check dimensions
if dim(obj) ~= dim(E)
    throw(CORAerror('CORA:dimensionMismatch',obj,E));
end

% simulatenous diagonalization: Find Tb such that
% Tb'*Q1*Tb = I and Tb'*Q2*Tb = D (diagonal)
% if max(diag(D))<=1 => contained
TOL = min(obj.TOL,E.TOL);
[~,D] = simdiag(obj.Q,E.Q,TOL);
tmp = max(diag(D));
res = tmp < 1+TOL | withinTol(tmp,1+TOL);

%------------- END OF CODE --------------