function res_empty = test_conZonotope_conZonotope
% test_conZonotope_conZonotope - unit test function of conZonotope (constructor)
%
% Syntax:  
%    res = test_conZonotope_conZonotope
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      19-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty conZonotope
conZono = conZonotope();
res_empty = true;
if ~isempty(conZono)
    res_empty = false;
end


if res_empty
    disp('test_conZonotope_conZonotope successful');
else
    disp('test_conZonotope_conZonotope failed');
end

%------------- END OF CODE --------------