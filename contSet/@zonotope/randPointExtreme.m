function [p] = randPointExtreme(obj)
% randPointExtreme - generates a random extreme points of a zonotope
%
% Syntax:  
%    [p] = randPointExtreme(obj)
%
% Inputs:
%    obj - zonotope object
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      14-May-2009
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%obtain number of generators
nrOfGenerators=length(obj.Z(1,:))-1;

%initialize the random point
p=obj.Z(:,1);

%add generators randomly
for i=1:nrOfGenerators
    val=sign(2*rand(1)-1);
    p=p+val*obj.Z(:,i+1);
end




%------------- END OF CODE --------------