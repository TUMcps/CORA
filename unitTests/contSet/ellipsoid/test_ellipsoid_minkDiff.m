function res = test_ellipsoid_minkDiff
% test_ellipsoid_minkDiff - unit test function of minkDiff
%
% Syntax:  
%    res = test_ellipsoid_minkDiff
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
    E1 = E_c{i}.E1;
    Ed1 = E_c{i}.Ed1;
    E0 = E_c{i}.E0;
    
    % since E0.Q = 0
    if isempty(minkDiff(E1,E0))
        res = false;
        break;
    end
    
    Ep = E0+Ed1;
    if ~contains(minkDiff(Ep,Ed1)+Ed1,Ep)
        res = false;
        break;
    end
    
end

%------------- END OF CODE --------------
