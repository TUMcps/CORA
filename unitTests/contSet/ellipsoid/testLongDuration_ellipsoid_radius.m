function res = testLongDuration_ellipsoid_radius
% testLongDuration_ellipsoid_radius - unit test function of radius
%
% Syntax:  
%    res = testLongDuration_ellipsoid_radius
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
% Written:      19-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty case: dim = 0
res_empty = true;
E = ellipsoid();
if rank(E) ~= 0
    res_empty = false;
end


% random tests
res_rand = true;
nrOfTests = 100;
bools = [false,true];
for i=1:nrOfTests
    
    % random dimension
    n = randi(30);
    for k=1:2
        E = ellipsoid.generateRandom(n,bools(k));
        E = ellipsoid(E.Q);
        [U,~,~] = svd(E.Q);
        IntE = interval(U'*E);
        R = rad(IntE);

        for j=1:dim(E)
            if ~all(withinTol(radius(E,j),R(1:j),E.TOL))
                res_rand = false;
                break;
            end
        end

        if ~res_rand
            break;
        end
    end
    if ~res_rand
        break;
    end
end

% combine results
res = res_empty && res_rand;

if res
    disp('testLongDuration_ellipsoid_radius successful');
else
    disp('testLongDuration_ellipsoid_radius failed');
end

%------------- END OF CODE --------------