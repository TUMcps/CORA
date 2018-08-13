function I = interval(Z)
% interval - Overapproximates a zonotope by an interval hull
%
% Syntax:  
%    I = interval(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    I - interval object
%
% Example: 
%    Z = zonotope(rand(2,5));
%    I = interval(Z);
%    plot(Z);
%    hold on
%    plot(I);
%
% Other m-files required: intervalhull(constructor)
% Subfunctions: none
% MAT-files required: none
%
% See also: vertices, polytope

% Author:       Matthias Althoff
% Written:      14-September-2006 
% Last update:  22-March-2007
%               11-November-2010
%               04-May-2011
%               22-July-2016
% Last revision: ---

%------------- BEGIN CODE --------------

%extract center
c = Z.Z(:,1);

%determine left and right limit
%specially designed for high performance
delta = sum(abs(Z.Z),2) - abs(c);
leftLimit = c - delta;
rightLimit = c + delta;

%instantiate interval hull
I = interval(leftLimit,rightLimit);

%------------- END OF CODE --------------