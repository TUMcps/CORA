function [p] = randPoint(obj)
% randPoint - generates a random point within a zonotope
%
% Syntax:  
%    [p] = randPoint(obj)
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

% Author: Matthias Althoff
% Written: 23-September-2008 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%obtain number of generators
nrOfGenerators=length(obj.Z(1,:))-1;

%initialize the random point
p=obj.Z(:,1);

%add generators randomly
for i=1:nrOfGenerators
    p=p+(2*rand(1)-1)*obj.Z(:,i+1);
end




%------------- END OF CODE --------------