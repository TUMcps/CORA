function res = containsPoint(Z,p,varargin)
% containsPoint - determines if the point p is inside the zonotope Z1
%
% Syntax:  
%    res = containsPoint(Z,p)
%    res = containsPoint(Z,p,tolerance)
%
% Inputs:
%    Z - zonotope object
%    p - point specified as a vector
%    tolerance - numerical tolerance up to which the point is allowed to 
%                outside the zonotope
%
% Outputs:
%    res - boolean whether the point is inside the zonotope or not
%
% Example: 
%    Z = zonotope([1;0],[1 0; 0 1]);
%    p = [1;0];
%    res = containsPoint(Z,p);
%    
%    plot(Z); hold on;
%    scatter(p(1),p(2),16,'r');
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

% parse input arguments
if nargin == 3
	tolerance = varargin{1};
else
	tolerance = 0;
end

% generate halfspace representation if empty
if isempty(Z.halfspace)
	Z = halfspace(Z);
end

%simple test: Is point inside the zonotope?
N = length(Z.halfspace.K);
inequality = (Z.halfspace.H*p - Z.halfspace.K <= tolerance * ones(N,1));

res = (all(inequality));

%------------- END OF CODE --------------
