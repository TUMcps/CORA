function res = testLongDuration_capsule_center
% testLongDuration_capsule_center - unit test function of center
%
% Syntax:  
%    res = testLongDuration_capsule_center
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
% Written:      28-August-2019
% Last update:  12-March-2021 (add empty case)
% Last revision:---

%------------- BEGIN CODE --------------


% random capsules
res = true;
for n=1:50
    
    % init capsule
    c_true = randn(n,1);
    C = capsule(c_true, randn(n,1), rand(1));

    % read center
    c = center(C);
    
    % check result
    if ~all(c == c_true)
        res = false; break;
    end
end

if res
    disp('test_center successful');
else
    disp('test_center failed');
end

%------------- END OF CODE --------------