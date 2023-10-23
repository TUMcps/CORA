function P = polytope(I)
% polytope - convert an interval to a polytope
%
% Syntax:
%    P = polytope(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    P - polytope object
%
% Example: 
%    I = interval([-1;0],[1;3]);
%    P = polytope(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/polytope

% Authors:       Niklas Kochdumper
% Written:       19-November-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check if the interval has the correct format
if size(I,2) > 1
    throw(CORAerror('CORA:wrongInputInConstructor',...
        'Input argument ''I'' has the wrong format!'));
end

% construct halfspace constraints
n = length(I);
C = [eye(n);-eye(n)];
d = [supremum(I);-infimum(I)];

% eliminate unbounded directions
idxInf = isinf(d);
C = C(~idxInf,:);
d = d(~idxInf);

% construct polytope object
P = polytope(C,d);

% ------------------------------ END OF CODE ------------------------------
