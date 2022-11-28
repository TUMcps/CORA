function Z = enclose(Z,varargin)
% enclose - encloses a zonotope and its affine transformation
%
% Description:
%    Computes the set
%    { a x1 + (1 - a) * (M x1 + x2) | x1 \in Z, x2 \in Z2, a \in [0,1] }
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

% Author:        Matthias Althoff
% Written:       30-September-2006 
% Last update:   22-March-2007
% Last revision: ---

%------------- BEGIN CODE --------------

% parse input arguments
if nargin == 2
    Z2 = varargin{1};
elseif nargin == 3
    Z2 = (varargin{1}*Z) + varargin{2};
else
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% retrieve number of generators of the zonotopes
generators1 = length(Z.Z(1,:));
generators2 = length(Z2.Z(1,:));

% if first zonotope has more or equal generators
if generators2 <= generators1
    Zcut = Z.Z(:,1:generators2);
    Zadd = Z.Z(:,(generators2+1):generators1);
    Zequal = Z2.Z;
else
    Zcut = Z2.Z(:,1:generators1);
    Zadd = Z2.Z(:,(generators1+1):generators2);
    Zequal = Z.Z;
end

% compute enclosing zonotope
Z.Z = [(Zcut+Zequal)/2,(Zcut-Zequal)/2,Zadd];

%------------- END OF CODE --------------