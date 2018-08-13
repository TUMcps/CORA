function [Z]=parallelotope(obj,Lambda)
% parallelotope - Computes a parallelotope in G-representation based on a
% coordinate transformation in which the transformed vertices are enclosed
% by an intervalhull
%
% Syntax:  
%    [Z]=parallelotope(obj,Lambda)
%
% Inputs:
%    obj - vertices object
%    Lambda - Lambda matrix determining the coordinate transformation
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      07-June-2009
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

%compute inverse of Lambda
invLambda=inv(Lambda);

%compute enclosing intervalhull in transformed coordinate system
IH=interval(invLambda*obj);
Z=Lambda*zonotope(IH);

%------------- END OF CODE --------------