function res = testLong_capsule_radius
% testLong_capsule_radius - unit test function of radius
%
% Syntax:
%    res = testLong_capsule_radius
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       12-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% random tests
res = true;
for n=1:2:30
    % init random capsule with generator of length zero
    r = rand(1);
    C = capsule(randn(n,1),zeros(n,1),r);
    
    % since capsule is a ball, enclosing radius is radius
    if ~withinTol(radius(C),r)
        res = false;
        break;
    end
end

% ------------------------------ END OF CODE ------------------------------
