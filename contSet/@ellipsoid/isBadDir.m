function res = isBadDir(E1,E2,L)
% isBadDir - checks if specified directions are bad directions for
%               Minkowski difference of E1 and E2
%
% Syntax:
%    res = isBadDir(E1,E2,L)
%
% Inputs:
%    L  - (n x N) matrix, where n must be dim(E1)=dim(E2), and N is the
%           number of directions to check
%    E1 - ellipsoid object
%    E2 - ellipsoid object
%
% Outputs:
%    res - boolean (vector; length(res)=size(L,2))
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       13-March-2019
% Last update:   10-June-2022
%                04-July-2022
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% input check
inputArgsCheck({{E1,'att','ellipsoid','scalar'};
                {E2,'att','ellipsoid','scalar'};
                {L,'att','numeric',{'nrows',dim(E1)}}});

% check dimension
if dim(E1)~=dim(E2)
    throw(CORAerror('CORA:wrongValue','second',...
                        '"E1" and "E2" need have the same dimension.'));
end

TOL = min(E1.TOL,E2.TOL);
[~,D] = simdiag(E1.Q,E2.Q,TOL);
r = 1/max(diag(D));
res = true(1,size(L,2));
for i=1:size(L,2)
    l = L(:,i);
    res(i) = sqrt(l'*E1.Q*l)/sqrt(l'*E2.Q*l) > r+TOL;
end

% ------------------------------ END OF CODE ------------------------------
