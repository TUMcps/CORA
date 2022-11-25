function intMat = randomIntervalMatrix(center,delta,varyingElements)
% randomIntervalMatrix - generates a random interval matrix with a
% specified center and a specified delta matrix or scalar. The number of
% elements of that matrix which are uncertain has to be specified, too
%
% Syntax:  
%     intMat = randomIntervalMatrix(center,delta,varyingElements)
%
% Inputs:
%    center - center matrix
%    delta - delta matrix
%    varyingElements - number of elements that may vary
%
% Outputs:
%    intMat - interval matrix
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
dim = length(center);

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
    
%add center
intMat = intervalMatrix(center,delta);

%------------- END OF CODE --------------