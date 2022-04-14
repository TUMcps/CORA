function [c] = center(E)
% center - returns the center of an ellipsoid object
%
% Syntax:  
%    [c] = center(E);
%
% Inputs:
%    E - Ellipsoid object
%
% Outputs:
%    c - center of E
%
% Example: 
%    E1=ellipsoid([1 0; 0 1]);
%    c = center(E1);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      13-March-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
c = E.q;
%------------- END OF CODE --------------