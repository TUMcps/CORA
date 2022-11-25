function res = test_halfspace_generateRandom
% test_halfspace_generateRandom - unit test function of generateRandom
%
% Syntax:  
%    res = test_halfspace_generateRandom
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
% Written:      27-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

try
    % instantiate halfspaces
    h = halfspace.generateRandom();
    % fix dimension
    h = halfspace.generateRandom(3);
    if dim(h) ~= 3
        error("Fixed dimension not obeyed");
    end
    res = true;
catch
    res = false;
end


if res
    disp('test_generateRandom successful');
else
    disp('test_generateRandom failed');
end

%------------- END OF CODE --------------