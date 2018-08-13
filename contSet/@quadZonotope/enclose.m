function [qZ] = enclose(qZ1,qZ2)
% enclose - Generates an overapproximative quadZonotope that encloses two 
% quadZonotopes of equal dimension
%
% Syntax:  
%    [qZ] = enclose(qZ1,qZ2)
%
% Inputs:
%    qZ1 - first quadZonotope object
%    qZ2 - second quadZonotope object
%
% Outputs:
%    qZ - quadZonotope enclosing qZ1 and qZ2
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      05-September-2012
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%overapproximate quadZonotopes by zonotopes
Z1 = zonotope(qZ1);
Z2 = zonotope(qZ2);

%compute enclosing zonotope
Z=enclose(Z1,Z2);

%write result as a quadZonotope
Zmat = get(Z,'Z');
c = Zmat(:,1);
if isempty(qZ1.G)
    len = 0;
    G = 0*c;
else
    len = length(qZ1.G(1,:));
    G = Zmat(:,2:len+1);
end
Grest = Zmat(:,len+2 : end);

qZ = quadZonotope(c,G,[],[],Grest);

%------------- END OF CODE --------------