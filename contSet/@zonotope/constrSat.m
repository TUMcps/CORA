function res = constrSat(obj,C,d)
% constrSat - checks if all values x of a zonotope satisfy the constraint
% Cx<=d
%
% Syntax:  
%    res = constrSat(obj,C,d)
%
% Inputs:
%    obj - zonotope object
%    C - normal vectors of constraints
%    d - distance to origin of constraints
%
% Outputs:
%    res - 1/0 if constraint is satisfied
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
IH = interval(C*obj) + (-d);

%check if interval contains 0
res = all(supremum(IH)<0);


%------------- END OF CODE --------------