function [Z1] = enclose(Z1,Z2)
% enclose - Generates a zonotope that encloses two zonotopes of equal order
% and dimension
%
% Syntax:  
%    [Z1] = enclose(Z1,Z2)
%
% Inputs:
%    Z1 - first zonotope object
%    Z2 - second zonotope object
%
% Outputs:
%    Z1 - zonotope, that encloses Z1 and Z2
%
% Example: 
%    Z1=zonotope(rand(2,4));
%    Z2=zonotope(rand(2,4));
%    Z=enclose(Z1,Z2);
%    plot(Z1);
%    hold on
%    plot(Z2);
%    plot(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 30-September-2006 
% Last update: 22-March-2007
% Last revision: ---

%------------- BEGIN CODE --------------

%retrieve number of generators of the zonotopes
generators1=length(Z1.Z(1,:));
generators2=length(Z2.Z(1,:));

%if first zonotope has more or equal generators
if generators2<=generators1
    Zcut=Z1.Z(:,1:generators2);
    Zadd=Z1.Z(:,(generators2+1):generators1);
    Zequal=Z2.Z;
else
    Zcut=Z2.Z(:,1:generators1);
    Zadd=Z2.Z(:,(generators1+1):generators2);
    Zequal=Z1.Z;
end

%compute enclosing zonotope
Z1.Z=[(Zcut+Zequal)/2,(Zcut-Zequal)/2,Zadd];

%------------- END OF CODE --------------