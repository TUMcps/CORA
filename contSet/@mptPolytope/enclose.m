function obj = enclose(varargin)
% enclose - Generates a polytope that encloses a polytope and its linear 
%           transformation
%
% Syntax:  
%    P = enclose(P1,P2)
%    P = enclose(P1,M,Pplus)
%
% Inputs:
%    P1 - first polytope object
%    P2 - second polytope object, satisfying P2 = (M * P1) + Pplus
%    M - matrix for the linear transformation
%    Pplus - polytope object added to the linear transformation
%
% Outputs:
%    obj - polytope that encloses P1 and P2
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclose

% Author:       Matthias Althoff
% Written:      02-February-2011
% Last update:  12-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

% parse input arguments
if nargin == 2
    obj = varargin{1};
    P = varargin{2};
else
    obj = varargin{1};
    M = varargin{2};
    Pplus = varargin{3};
    
    P = (M*obj) + Pplus;
end

% call MPT-toolbox method to over-approximate with convex hull
try %MPT3
    obj.P = convexHull(obj.P, P.P);
catch %MPT2
    %compute convex hull
    obj.P = hull([obj.P P.P]);
end


%------------- END OF CODE --------------