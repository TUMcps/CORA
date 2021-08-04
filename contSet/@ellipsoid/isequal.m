function res = isequal(E1,E2)
% isequal - checks if ellipsoids E1, E2 are equal within tolerance
%
% Syntax:  
%    B = isequal(E,Y) gives a logical indicating whether E1,E2 are equal
%    (within tolerance)
%
% Inputs:
%    E1 - ellipsoids object
%    E2 - ellipsoids object
%
% Outputs:
%    res - boolean value indicating whether points E1,E2 are equal
%
% Example: 
%    E1 = ellipsoid([1,0;0,1/2],[1;1]);
%    E2 = ellipsoid([1+1e-15,0;0,1/2],[1;1]);
%    B = equal(E1,E2);
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
% Last revision:---

%------------- BEGIN CODE --------------
if ~isa(E1,'ellipsoid') || ~isa(E2,'ellipsoid')
    error('Both input arguments need to be of type "ellipsoid"!');
end
res = E1==E2;
%------------- END OF CODE --------------