function [M] = mapFromDirection(direction)
% mapFromDirection - computes a linear map for parallelotope enclosure 
% given the direction of the flow field
%
% Syntax:  
%    [M] = mapFromDirection(direction)
%
% Inputs:
%    direction - vector specifying the direction
%
% Outputs:
%    M - matrix specifying the linear map
%
% Example: 
%    ---
%
% Other m-files required: ---
% Subfunctions: none
% MAT-files required: none
%
% See also: interval,  vertices

% Author:       Matthias Althoff
% Written:      26-October-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

dim=length(direction);
orient=eye(dim);

%check if direction is not the origin
if norm(direction)>0
    newGen=direction/norm(direction);
else
    newGen=[1;zeros(dim-1,1)];
end

%retrieve most aligned generator from orient
for iGen=1:length(orient(1,:))
     h(iGen)=abs(newGen'*orient(:,iGen)/norm(orient(:,iGen)));
end

[val,ind]=sort(h);
pickedIndices=ind(1:(end-1));

M=[newGen,orient(:,pickedIndices)];
 
%------------- END OF CODE --------------