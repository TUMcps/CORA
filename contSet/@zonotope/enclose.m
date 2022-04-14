function Z = enclose(varargin)
% enclose - returns a zonotope enclosing a zonotope and its linear 
%    transformation
%
% Syntax:  
%    Z = enclose(Z1,Z2)
%    Z = enclose(Z1,M,Zplus)
%
% Inputs:
%    Z1 - first zonotope object
%    Z2 - second zonotope object, satisfying Z2 = (M * Z1) + Zplus
%    M - matrix for the linear transformation
%    Zplus - zonotope object added to the linear transformation
%
% Outputs:
%    Z - zonotope, that encloses Z1 and Z2
%
% Example: 
%    Z1 = zonotope([1.5 1 0; 1.5 0 1]);
%    Z2 = [-1 0; 0 -1]*Z1 + [0.5;0.5];
%
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
    Z1 = varargin{1};
    Z2 = varargin{2};
else
    Z1 = varargin{1};
    M = varargin{2};
    Zplus = varargin{3};
    
    Z2 = (M*Z1) + Zplus;
end

% retrieve number of generators of the zonotopes
generators1=length(Z1.Z(1,:));
generators2=length(Z2.Z(1,:));

% if first zonotope has more or equal generators
if generators2<=generators1
    Zcut=Z1.Z(:,1:generators2);
    Zadd=Z1.Z(:,(generators2+1):generators1);
    Zequal=Z2.Z;
else
    Zcut=Z2.Z(:,1:generators1);
    Zadd=Z2.Z(:,(generators1+1):generators2);
    Zequal=Z1.Z;
end

% compute enclosing zonotope
Z = zonotope([(Zcut+Zequal)/2,(Zcut-Zequal)/2,Zadd]);

%------------- END OF CODE --------------