function res = test_capsule_generateRandom
% test_capsule_generateRandom - unit test function of generateRandom
%
% Syntax:  
%    res = test_capsule_generateRandom
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

% run through all possibilities
try
    % random dim and cen
    C = capsule.generateRandom();
    % fixed dim, random cen
    C = capsule.generateRandom(3);
    if dim(C) ~= 3
        error("Fixed dimension not obeyed");
    end
    % fixed dim and cen
    C = capsule.generateRandom(2,[1;0]);
    if dim(C) ~= 2 && ~isequal(center(C),[1;0])
        error("Fixed dimension and center not obeyed");
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