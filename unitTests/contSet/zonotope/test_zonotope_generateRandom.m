function res = test_zonotope_generateRandom
% test_zonotope_generateRandom - unit test function of generateRandom
%
% Syntax:  
%    res = test_zonotope_generateRandom
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
    Z = zonotope.generateRandom();
    % fixed dim, random cen
    Z = zonotope.generateRandom(3);
    if dim(Z) ~= 3
        error("Fixed dimension not obeyed");
    end
    % fixed dim and cen
    Z = zonotope.generateRandom(2,[1;0]);
    if dim(Z) ~= 2 && ~isequal(center(Z),[1;0])
        error("Fixed dimension and center not obeyed");
    end
    % fixed dim, cen and gen
    Z = zonotope.generateRandom(3,[1;0;-3],35);
    if dim(Z) ~= 3 && ~isequal(center(Z),[1;0;-3]) && ...
            length(generators(Z)) == 35
        error("Fixed dimension, center and number of generators not obeyed");
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