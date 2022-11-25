function [Zbundle]=reduce(Zbundle,option,varargin)
% reduce - Reduces the order of a zonotope bundle; see reduce zonotope
%
% Syntax:  
%    [ZbundleRed,t]=reduce(Zbundle,option,order,filterLength)
%
% Inputs:
%    Zbundle - zonotope bundle
%    option - reduction method selector
%    order - maximum order of reduced zonotope
%
% Outputs:
%    Zbundle - bundle of reduced zonotopes
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
% Written:      09-November-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%reduce for each zonotope
for i=1:Zbundle.parallelSets
    Zbundle.Z{i}=reduce(Zbundle.Z{i},option,varargin{:});
end

%------------- END OF CODE --------------