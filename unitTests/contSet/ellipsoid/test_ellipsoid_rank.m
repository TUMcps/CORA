function res = test_ellipsoid_rank
% test_ellipsoid_rank - unit test function of rank
%
% Syntax:  
%    res = test_ellipsoid_rank
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
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
load cases.mat E_c
for i=1:length(E_c)
    E1 = E_c{i}.E1; % non-deg
    Ed1 = E_c{i}.Ed1; % deg
    E0 = E_c{i}.E0; % all zero
    n = dim(E1);
    
    if rank(E1)~=n || rank(Ed1)==n || rank(Ed1)<=rank(E0) || rank(E0)~=0
        res = false;
        break;
    end
    
end


if res
    disp([mfilename,' successful']);
else
    disp([mfilename,' failed']);
end
%------------- END OF CODE --------------
