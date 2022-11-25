function res = isFullDim(obj)
% isFullDim - check if a constrained zonotope is full-dimensional
%
% Syntax:  
%    res = isFullDim(obj)
%
% Inputs:
%    obj - conZonotope object
%
% Outputs:
%    res - 1 if conZonotope is full-dimensional, 0 else
%
% Example:
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1];
%    b = 1;
%    cZono1 = conZonotope(Z,A,b);
%
%    hp = conHyperplane([1,-2],1);
%    cZono2 = cZono1 & hp;
%
%    isFullDim(cZono1)
%    isFullDim(cZono2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isFullDim

% Author:       Niklas Kochdumper
% Written:      02-January-2020 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    if isempty(obj.A)
        
        % call zonotope isFullDim method
        res = isFullDim(zonotope(obj.Z));
        
    else

        % compute null-space of the constraints
        T = null(obj.A);

        % transform generator matrix into the null-space
        G_ = obj.Z(:,2:end) * T;

        % check if rank of generator matrix is equal to the dimension
        dimG = size(G_,1);
        rankG = rank(G_);

        res = dimG == rankG;    
    end

%------------- END OF CODE --------------