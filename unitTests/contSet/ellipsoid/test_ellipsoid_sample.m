function res = test_ellipsoid_sample
% test_ellipsoid_sample - unit test function of in
%
% Syntax:  
%    res = test_ellipsoid_sample
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

% Author:       Victor Gaßmann
% Written:      15-October-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
nRuns = 2;
%supress warning if ellipsoid is too squished (in sample(E,1000))
warning('off','all');
for i=10:5:20
    for j=1:nRuns
        E = ellipsoid.generateRandom(i);
        Y = sample(E,100);
        for k=1:length(Y)
            if (Y(:,k)-E.q)'*inv(E.Q)*(Y(:,k)-E.q)>(1+E.TOL)
                res = false;
                break;
            end
        end
    end
    if ~res
        break;
    end
end
warning('on','all');
if res
    disp('test_ellipsoid_sample successful');
else
    disp('test_ellipsoid_sample failed');
end
%------------- END OF CODE --------------
