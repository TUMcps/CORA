function result = containsPoint(Z1,p)
% containsPoint - determines if the point p is inside the zonotope Z1
%
% Syntax:  
%    result = containsPoint(Z1,p)
%
% Inputs:
%    Z1 - zonotope object
%    p - point specified as a vector
%
% Outputs:
%    result - 1/0 if point is inside the zonotope or not
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Niklas Kochdumper
% Written:      30-January-2018 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% generate halfspace representation if empty
if isempty(Z1.halfspace)
    Z1 = halfspace(Z1);
end

%simple test: Is point inside the zonotope?
inequality = (Z1.halfspace.H*p<=Z1.halfspace.K);

result = (all(inequality));

%------------- END OF CODE --------------
