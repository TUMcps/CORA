function [Zbundle]=and(Zbundle1,Zbundle2)
% and - returns the intersection of two zonotope bundles or a zonotope
% bundle and a zonotope
%
% Syntax: 
%    [Zbundle1]=and(Zbundle1,Zbundle2)
%
% Inputs:
%    Zbundle1 - 1st zonotope (bundle)
%    Zbundle1 - 2nd zonotope (bundle)
%
% Outputs:
%    Zbundle - zonotope bundle after intersection
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: ---
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      16-November-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%Find a zonotope bundle object
%Is summand1 a zonotope bundle?
if isa(Zbundle1,'zonotopeBundle')
    %initialize resulting zonotope
    Zbundle=Zbundle1;
    %initialize other summand
    otherObj=Zbundle2;
%Is summand2 a zonotope bundle?    
elseif isa(Zbundle2,'zonotopeBundle')
    %initialize resulting zonotope
    Zbundle=Zbundle2;
    %initialize other summand
    otherObj=Zbundle1;
end

%check if other object is a zonotope or zonotope bundle
if isa(otherObj,'zonotope')
    Zbundle.Z{end+1} = otherObj;
    Zbundle.parallelSets = Zbundle.parallelSets+1;
elseif isa(otherObj,'zonotopeBundle')
    for i=1:otherObj.parallelSets
        Zbundle.Z{end+1} = otherObj.Z{i};
    end
    Zbundle.parallelSets = Zbundle.parallelSets + otherObj.parallelSets;
end


%------------- END OF CODE --------------