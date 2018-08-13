function [Zbundle1] = enclose(Zbundle1,Zbundle2)
% enclose - Generates a zonotope bundle that encloses two zonotopes bundles
% pairwise such that each bundle has to have equal number of zonotopes
%
% Syntax:  
%    [Zbundle1] = enclose(Zbundle1,Zbundle2)
%
% Inputs:
%    Zbundle1 - first zonotope bundle
%    Zbundle2 - second zonotope bundle
%
% Outputs:
%    Zbundle1 - zonotope bundle, that encloses both bundles
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
% Written:      10-November-2010 
% Last update:  25-January-2016
% Last revision:---

%------------- BEGIN CODE --------------

if isa(Zbundle2,'zonotopeBundle')
    %compute enclosure for each zonotope pair
    for i=1:Zbundle1.parallelSets
        Zbundle1.Z{i}=enclose(Zbundle1.Z{i},Zbundle2.Z{i});
    end
elseif isa(Zbundle2,'zonotope')
    %compute enclosure for each zonotope
    for i=1:Zbundle1.parallelSets
        Zbundle1.Z{i}=enclose(Zbundle1.Z{i},Zbundle2);
    end
end


%------------- END OF CODE --------------