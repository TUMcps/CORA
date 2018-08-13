function [res] = inORH(Z1,Z2,rotMatrixInv)
% inORH - checks if a zonotope Z1 is in an oriented rectangular hull (Z2)
% which has a certain rotation matrix
%
% Syntax:  
%    [res] = inORH(Z1,Z2,rotMatrix)
%
% Inputs:
%    Z1 - zonotope object
%    Z2 - ORH in zonotope representation
%    rotMatrixInv - inverse rotation matrix of the ORH
%
% Outputs:
%    res - 1/0 depending if Z1 is enclosed in Z2 or not
%
% Example: 
%
% Other m-files required: 
% Subfunctions: none
% MAT-files required: none
%
% See also: in

% Author: Matthias Althoff
% Written: 15-January-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%rotate Z1 and Z2
Z1rot=rotMatrixInv*Z1;
Z2rot=rotMatrixInv*Z2;

%convert rotated zonotopes to interval hulls
IH1=interval(Z1rot);
IH2=interval(Z2rot);

%check if IH1 in IH2
res=(IH1<=IH2);
    
%------------- END OF CODE --------------