function [cZ] = enclose(cZ1,cZ2)
% enclose - generates a conZonotope object that encloses two constrained 
%           zonotopes, where the second constrained zonotope is a linear 
%           transformation of the first one
%
% Syntax:  
%    [cZ] = enclose(cZ1, cZ2)
%
% Description: 
%    ATTENTION: This function can only be applied if the second conZonotope
%    object cZ2 is a linear transformation of the first conZonotope object
%    cZ1: 
%           cZ2 = T*cZ1 + t
%    where T is a numerical matrix and is a vector, a zonotope object or
%    an interval object
%
% Inputs:
%    cZ1 - first conZonotope object
%    cZ2 - second conZonotope object
%
% Outputs:
%    cZ - conZonotope object that encloses cZ1 and cZ2
%
% Example: 
%    % constrained zonotope
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1];
%    b = 1;
%    cZono1 = conZonotope(Z,A,b);
%
%    % linear transformation
%    T = [1 1.5;-0.5 1];
%    t = [4;7];
%    cZono2 = T*cZono1 + t;
%
%    % convex hull
%    cZonoRes = enclose(cZono1,cZono2);
%
%    % visualization
%    hold on
%    plotFilled(cZonoRes,[1,2],[0.6,0.6,0.6]);
%    plotFilled(cZono1,[1,2],'r');
%    plotFilled(cZono2,[1,2],'b');
%   
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Niklas Kochdumper
% Written: 28-June-2018 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

% retrieve number of generators of the zonotopes
g1 = length(cZ1.Z(1,:));
g2 = length(cZ2.Z(1,:));

% check if the constraints of both zonotopes are identical (otherwise it is
% not possible that they are linear transformations of each other)
if ~isempty(cZ1.A)
    
    g = min(g1,g2)-1;
    A1 = cZ1.A(:,1:g);
    A2 = cZ2.A(:,1:g);
    
    if ~all(size(A1) == size(A2)) || max(max(abs(A1-A2))) > eps || max(abs(cZ1.b-cZ2.b)) > eps
        error('Operation only vailid if conZonotope two is a linear transformation of conZonotope one!');
    end
elseif ~isempty(cZ2.A)
    error('Operation only vailid if conZonotope two is a linear transformation of conZonotope one!');
end

% divide generators into blocks
if g2 <= g1
    Z1 = cZ1.Z(:,1:g2);
    Zadd = cZ1.Z(:,(g2+1):end);
    Z2 = cZ2.Z;
    A = cZ2.A;
else
    Z2 = cZ2.Z(:,1:g1);
    Zadd = cZ2.Z(:,(g1+1):end);
    Z1 = cZ1.Z;
    A = cZ1.A;
end

% construct enclosing constrained zonotope
cZ = cZ1;

cZ.Z = [(Z1+Z2)/2, (Z1-Z2)/2, Zadd];
cZ.A = [A, zeros(size(A,1),size(Z1,2) + size(Zadd,2))];

cZ.R = [];
cZ.ksi = [];

%------------- END OF CODE --------------