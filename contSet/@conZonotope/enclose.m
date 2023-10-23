function cZ = enclose(cZ,varargin)
% enclose - encloses a constrained zonotope and its affine transformation
%
% Description:
%    Computes the set
%    { a x1 + (1 - a) * x2 | x1 \in cZ, x2 \in cZ2, a \in [0,1] }
%    where cZ2 = M*cZ + cZplus
%
% Syntax:
%    cZ = enclose(cZ,cZ2)
%    cZ = enclose(cZ,M,cZplus)
%
% Inputs:
%    cZ - conZonotope object
%    cZ2 - conZonotope object
%    M - matrix for the linear transformation
%    cZplus - conZonotope object added to the linear transformation
%
% Outputs:
%    cZ - conZonotope object that encloses cZ and cZ2
%
% Example: 
%    % constrained zonotope
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1]; b = 1;
%    cZ1 = conZonotope(Z,A,b);
%
%    % linear transformation
%    M = [1 1.5;-0.5 1];
%    cZplus = [4;7];
%    cZ2 = M*cZ1 + cZplus;
%
%    % enclose
%    cZonoRes = enclose(cZ1,cZ2);
%
%    % visualization
%    figure; hold on;
%    plot(cZonoRes,[1,2],'FaceColor',[0.6,0.6,0.6]);
%    plot(cZ1,[1,2],'FaceColor','r');
%    plot(cZ2,[1,2],'FaceColor','b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclose

% Authors:       Niklas Kochdumper
% Written:       28-June-2018 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
if nargin == 2
    cZ2 = varargin{1};
elseif nargin == 3
    % compute M*cZ1 + cZplus
    cZ2 = (varargin{1}*cZ) + varargin{2};
else
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% retrieve number of generators of the zonotopes
g1 = size(cZ.G,2);
g2 = size(cZ2.G,2);

% check if the constraints of both zonotopes are identical (otherwise it is
% not possible that they are linear transformations of each other)
if ~isempty(cZ.A)
    
    g = min(g1,g2);
    A1 = cZ.A(:,1:g);
    A2 = cZ2.A(:,1:g);
    
    if ~all(size(A1) == size(A2)) || max(max(abs(A1-A2))) > eps || max(abs(cZ.b-cZ2.b)) > eps
        throw(CORAerror('CORA:specialError',...
            'Operation only valid if conZonotope two is a linear transformation of conZonotope one!'));
    end
elseif ~isempty(cZ2.A)
    throw(CORAerror('CORA:specialError',...
        'Operation only valid if conZonotope two is a linear transformation of conZonotope one!'));
end

% divide generators into blocks
if g2 <= g1
    cG = (cZ.c-cZ2.c)/2;
    G1 = cZ.G(:,1:g2);
    Gadd = cZ.G(:,(g2+1):end);
    G2 = cZ2.G;
    A = cZ2.A;
else
    cG = (cZ2.c-cZ.c)/2;
    G2 = cZ2.G(:,1:g1);
    Gadd = cZ2.G(:,(g1+1):end);
    G1 = cZ.G;
    A = cZ.A;
end

% construct enclosing constrained zonotope
cZ.c = (cZ.c+cZ2.c)/2;
cZ.G = [(G1+G2)/2, cG, (G1-G2)/2, Gadd];
cZ.A = [A, zeros(size(A,1),size(G1,2) + 1 + size(Gadd,2))];

cZ.R = [];
cZ.ksi = [];

% ------------------------------ END OF CODE ------------------------------
