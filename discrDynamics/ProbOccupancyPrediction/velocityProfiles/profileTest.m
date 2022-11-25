function profileTest(handle,x0,xEnd,acc)
% profileTest - tests velocity profiles
%
% Syntax:  
%    profileTest(handle,x0,xEnd)
%
% Inputs:
%    handle - handle of the profile function
%    x0 - beginning position
%    xEnd - end position
%
% Outputs:
%    ---
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 03-July-2008 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%generate position vector
posVector=x0:1:xEnd;

%obtain velocities
for i=1:length(posVector)
    velVector(i)=handle(posVector(i),acc);
end

%plot result
plot(posVector,velVector);

%------------- END OF CODE --------------