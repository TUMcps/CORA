function E = ellipsoid(obj)
% ellipsoid - Converts an interval object into an ellipsoid object
%
% Syntax:  
%    Z = ellipsoid(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    I = interval([1;-1], [2; 1]);
%    E = ellipsoid(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Author:       Victor Gaﬂmann
% Written:      15-October-2019 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
%convert interval to zonotope
Z = zonotope(obj);
%n==m, therefore exact computation very efficient
E = ellipsoid(Z,'o:exact');
%------------- END OF CODE --------------