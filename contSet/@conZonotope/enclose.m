function cZ = enclose(varargin)
% enclose - generates a conZonotope object that encloses two constrained 
%           zonotopes, where the second constrained zonotope is a linear 
%           transformation of the first one
%
% Syntax:  
%    cZ = enclose(cZ1, cZ2)
%    cZ = enclose(cZ1, M, cZplus)
%
% Inputs:
%    cZ1 - first conZonotope object
%    cZ2 - second conZonotope object, satisfying cZ2 = (M * cZ1) + cZplus
%    M - matrix for the linear transformation
%    cZplus - conZonotope object added to the linear transformation
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
%    M = [1 1.5;-0.5 1];
%    cZplus = [4;7];
%    cZono2 = M*cZono1 + cZplus;
%
%    % convex hull
%    cZonoRes = enclose(cZono1,cZono2);
%
%    % visualization
%    hold on
%    plot(cZonoRes,[1,2],'FaceColor',[0.6,0.6,0.6],'Filled',true);
%    plot(cZono1,[1,2],'r','Filled',true);
%    plot(cZono2,[1,2],'b','Filled',true);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclose

% Author: Niklas Kochdumper
% Written: 28-June-2018 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

% parse input arguments
if nargin == 2
    cZ1 = varargin{1};
    cZ2 = varargin{2};
else
    cZ1 = varargin{1};
    M = varargin{2};
    cZplus = varargin{3};
    
    cZ2 = (M*cZ1) + cZplus;
end

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