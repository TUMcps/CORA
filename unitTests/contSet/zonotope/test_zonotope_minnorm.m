function res = test_zonotope_minnorm
% test_zonotope_minnorm - unit test function of minnorm
%
% Syntax:  
%    res = test_zonotope_minnorm
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
% Written:      15-October-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
TOL = 1e-12;
res = true;
for i=2:7
    for j=i:5:20
        Z = zonotope([zeros(i,1),randn(i,j)]);
        val = minnorm(Z);
        %check for random unit directions if val<=suppfnc(Z,l)
        L = eq_point_set(i-1,2*j);
        for k=1:length(L)
            if (val-supportFunc(Z,L(:,k)))>TOL
                res = false;
                break;
            end
        end
        if ~res
            break;
        end
    end
    if ~res
        break;
    end
end
if res
    disp('test_zonotope_minnorm successful');
else
    disp('test_zonotope_minnorm failed');
end

%------------- END OF CODE --------------
