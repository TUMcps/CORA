function [IH] = interval(Z)
% interval - Overapproximates a zonotope by an interval
%
% Syntax:  
%    [IH] = interval(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    IH - interval object
%
% Example: 
%    Z=zonotope(rand(2,5));
%    IH=interval(Z);
%    plot(Z);
%    hold on
%    plot(IH);
%
% Other m-files required: intervalhull(constructor)
% Subfunctions: none
% MAT-files required: none
%
% See also: vertices, polytope

% Author:       Matthias Althoff
% Written:      14-September-2006 
% Last update:  22-March-2007
%               25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%extract generators
G=Z.Z(:,2:end);

%extract center
c=Z.Z(:,1);

%determine left and right limit
leftLimit=c-sum(abs(G),2);
rightLimit=c+sum(abs(G),2);

%instantiate interval hull
IH=interval(leftLimit,rightLimit);

%------------- END OF CODE --------------