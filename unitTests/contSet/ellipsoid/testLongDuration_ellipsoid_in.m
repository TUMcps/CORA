function res = testLongDuration_ellipsoid_in
% testLongDuration_ellipsoid_in - unit test function of in
%
% Syntax:  
%    res = testLongDuration_ellipsoid_in
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
% Last update:  07-August-2020
%               19-March-2021
% Last revision:---

%------------- BEGIN CODE --------------
res = true;

% check ellipsoid 
res = res && testLongDuration_component_ellipsoid_inEllipsoid;

% check zonotope
res = res && testLongDuration_component_ellipsoid_inZonotope;

% check point containment
runs = 10;
res_point = true;
bools = [false,true];
for i=1:5:10
    for j=1:runs
        N = 10*i;
        for k=1:2
            E = ellipsoid.generateRandom(i,bools(k));
            samples = randPoint(E,N);
            if ~in(E,samples)
                res_point = false;
                break;
            end
        end
        if ~res_point
            break;
        end
    end
    if ~res_point
        break;
    end
end

res = res && res_point;

if res
    disp('testLongDuration_ellipsoid_in successful');
else
    disp('testLongDuration_ellipsoid_in failed');
end
%------------- END OF CODE --------------