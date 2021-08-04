function res = testLongDuration_zonotope_supportFunc
% testLongDuration_zonotope_supportFunc - unit test function of support function
%
% Syntax:  
%    res = testLongDuration_zonotope_supportFunc
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

% Author:       Victor Gassmann
% Written:      11-October-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
dims = 2:5;
dGen = 5;
steps = 3;
res = true;
TOL = 1e-8;
for i=dims
    n = i;
    for j=1:steps
        m = n + j*dGen;
        %randomly generate zonotope 
        Z = zonotope([randn(n,1),randn(n,m)]);
        %check if [~,x]=minnorm(Z)^2==supportFunc(Z,x)
        %check if [~,x]=norm(Z,2,'exact')^2==supportFunc(Z,x)
        [val_min,x_min] = minnorm(Z);
        [val_max,x_max] = norm(Z,2,'exact');
        if abs(supportFunc(Z,x_min)-val_min^2)>TOL || ...
                abs(supportFunc(Z,x_max)-val_max^2)>TOL
            res = true;
            break;
        end
    end
    if ~res
        break;
    end
end
if res
    disp('testLongDuration_zonotope_supportFunc successful');
else
    disp('testLongDuration_zonotope_supportFunc failed');
end

%------------- END OF CODE --------------