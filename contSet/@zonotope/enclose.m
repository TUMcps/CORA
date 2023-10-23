function Z = enclose(Z,varargin)
% enclose - encloses a zonotope and its affine transformation
%
% Description:
%    Computes the set
%    { a x1 + (1 - a) * x2 | x1 \in Z, x2 \in Z2, a \in [0,1] }
%    where Z2 = M*Z + Zplus
%
% Syntax:
%    Z = enclose(Z,Z2)
%    Z = enclose(Z,M,Zplus)
%
% Inputs:
%    Z - zonotope object
%    Z2 - zonotope object
%    M - matrix for the linear transformation
%    Zplus - zonotope object added to the linear transformation
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    Z1 = zonotope([1.5 1 0; 1.5 0 1]);
%    M = [-1 0; 0 -1];
%    Z2 = M*Z1 + [0.5;0.5];
%    Z=enclose(Z1,Z2);
%
%    figure; hold on;
%    plot(Z1,[1,2],'r');
%    plot(Z2,[1,2],'g');
%    plot(Z,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/enclose

% Authors:       Matthias Althoff
% Written:       30-September-2006 
% Last update:   22-March-2007
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
if nargin == 2
    Z2 = varargin{1};
elseif nargin == 3
    M = varargin{1};
    Zplus = varargin{2};
    Z2 = (M*Z) + Zplus;
else
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% retrieve number of generators of the zonotopes
generators1 = size(Z.G,2);
generators2 = size(Z2.G,2);

% if first zonotope has more or equal generators
if generators2 <= generators1
    cG = (Z.c-Z2.c)/2;
    Gcut = Z.G(:,1:generators2);
    Gadd = Z.G(:,(generators2+1):generators1);
    Gequal = Z2.G;
else
    cG = (Z2.c-Z.c)/2;
    Gcut = Z2.G(:,1:generators1);
    Gadd = Z2.G(:,(generators1+1):generators2);
    Gequal = Z.G;
end

% compute enclosing zonotope
Z.c = (Z.c+Z2.c)/2;
Z.G = [(Gcut+Gequal)/2,cG,(Gcut-Gequal)/2,Gadd];

% ------------------------------ END OF CODE ------------------------------
