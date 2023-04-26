function deltaRand = randomDelta(delta,varyingElements)
% randomDelta - generates a random delta matrix with specified bounds; the 
%    number of uncertain elements of that matrix can also be specified
%
% Syntax:  
%     deltaRand = randomDelta(delta,varyingElements)
%
% Inputs:
%    delta - delta matrix
%    varyingElements - number of elements that may vary
%
% Outputs:
%    deltaRand - random delta matrix
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      25-June-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%determine dimension
dim = length(delta);

%check if varying elements are proper set
if varyingElements>dim^2
    varyingElements = dim^2;
end

%initialize random delta matrix
deltaRand = zeros(dim);

%obtain uncertain entries of delta
counter = 0;
while counter<varyingElements
    %get random position of varying element
    n = randi(dim);
    m = randi(dim);
    
    %check if element is already non-zero
    if deltaRand(n,m)==0
        deltaRand(n,m) = rand(1)*delta(n,m);
        counter = counter+1;
    end
end

%------------- END OF CODE --------------