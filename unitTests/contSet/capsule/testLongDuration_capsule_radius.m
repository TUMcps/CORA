function res = testLongDuration_capsule_radius
% testLongDuration_capsule_radius - unit test function of radius
%
% Syntax:  
%    res = testLongDuration_capsule_radius
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
% Written:      12-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% 2. Random tests
res = true;
for n=1:2:30
    % init random capsule with generator of length zero
    r = rand(1);
    C = capsule(randn(n,1),zeros(n,1),r);
    
    % since capsule is a ball, enclosing radius is radius
    if radius(C) ~= r
        res = false; break;
    end
end


% combine results
if res
    disp('test_radius successful');
else
    disp('test_radius failed');
end

%------------- END OF CODE --------------