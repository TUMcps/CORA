function res = test_ellipsoid_containsPoint
% test_ellipsoid_containsPoint - unit test function of containsPoint
%
% Syntax:  
%    res = test_ellipsoid_containsPoint
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

% Author:       Victor Gaﬂmann, Matthias Althoff
% Written:      13-March-2019
% Last update:  02-Sep-2019 (rename containsPoint)
% Last revision:---

%------------- BEGIN CODE --------------

%NOTICE: Before executing this test, make sure that test_ellipsoid_supportFunc
%works since this test uses supportFunc.

N = 1000;
E = ellipsoid.generateRandom();
samples = sample(E,N);
if ~all(containsPoint(E,samples))
    disp('test_ellipsoid_containsPoint failed');
    res = false;
else
    disp('test_ellipsoid_containsPoint successful');
    res = true;
end

%------------- END OF CODE --------------