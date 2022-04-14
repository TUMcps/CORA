function res = test_interval_generateRandom
% test_interval_generateRandom - unit test function of generateRandom
%
% Syntax:  
%    res = test_interval_generateRandom
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
    I = interval.generateRandom();
    % fixed dim, random cen
    I = interval.generateRandom(3);
    if dim(I) ~= 3
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