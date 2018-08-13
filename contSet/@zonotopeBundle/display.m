function display(obj)
% display - Displays a zonotope bundle
%
% Syntax:  
%    display(obj)
%
% Inputs:
%    Z - zonotope bundle
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

% Author:       Matthias Althoff
% Written:      09-November-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%display each zonotope
for i=1:obj.parallelSets
    disp(['zonotope ',num2str(i),':']);
    display(obj.Z{i});
end

%------------- END OF CODE --------------