function res = constrSat(Z,C,d)
% constrSat - checks if all points x within a zonotope satisfy the
%    constraint Cx <= d
%
% Syntax:  
%    res = constrSat(Z,C,d)
%
% Inputs:
%    Z - zonotope object
%    C - normal vectors of constraints
%    d - distance to origin of constraints
%
% Outputs:
%    res - boolean whether constraint is satisfied
%
% Example: 
%
% Other m-files required: 
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      10-August-2011
% Last update:  14-May-2017
% Last revision:---

%------------- BEGIN CODE --------------

%check if constraints are violated
IH = interval(C*Z) + (-d);

%check if interval contains 0
res = all(supremum(IH)<0);


%------------- END OF CODE --------------