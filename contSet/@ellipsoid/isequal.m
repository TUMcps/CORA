function res = isequal(E,S)
% isequal - checks if two ellipsoids are equal
%
% Syntax:  
%    res = isequal(E,S)
%
% Inputs:
%    E - ellipsoid object
%    S - ellipsoid object
%
% Outputs:
%    res - true/false
%
% Example: 
%    E1 = ellipsoid([1,0;0,1/2],[1;1]);
%    E2 = ellipsoid([1+1e-15,0;0,1/2],[1;1]);
%    res = isequal(E1,E2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      13-March-2019
% Last update:  15-October-2019
%               19-March-2021 (use 'eq')
%               04-July-2022 (VG: class array case)
% Last revision:---

%------------- BEGIN CODE --------------

inputArgsCheck({{E,'att','ellipsoid','scalar'};
                {S,'att','ellipsoid','scalar'}});

res = E==S;

%------------- END OF CODE --------------