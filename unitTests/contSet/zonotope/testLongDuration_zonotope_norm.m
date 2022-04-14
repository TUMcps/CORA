function res = testLongDuration_zonotope_norm
% testLongDuration_zonotope_norm - unit test function of norm
%
% Syntax:  
%    res = testLongDuration_zonotope_norm
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
% Written:      31-July-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
TOL = 1e-6;
res = true;
for i=2:4
    for j=i:2:10
        Z = zonotope([zeros(i,1),10*randn(i,j)]);
        % 2 norm test
        val2_exact = norm(Z,2,'exact');
        val2_ub = norm(Z,2,'ub');
        val2_ubc = norm(Z,2,'ub_convex');
        V = vertices(Z);
        if (val2_exact-val2_ub)>TOL || (val2_exact-val2_ubc)>TOL ...
                || abs(val2_exact-max(sqrt(sum(V.^2))))/val2_exact>TOL
            res = false;
        end
        % rest case are simply norm of interval
    end
    if ~res
        break;
    end
end
if res
    disp('testLongDuration_zonotope_norm successful');
else
    disp('testLongDuration_zonotope_norm failed');
end

%------------- END OF CODE --------------
