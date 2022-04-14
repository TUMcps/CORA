function res = testLongDuration_ellipsoid_cartProd
% testLongDuration_ellipsoid_cartProd - unit test function of cartProd
%
% Syntax:  
%    res = testLongDuration_ellipsoid_cartProd
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
% Written:      18-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

runs = 5;
res = true;
bools = logical([1,1,0,0;1,0,1,0]);
for j=1:runs
    for k=1:4
        E1 = ellipsoid.generateRandom(bools(1,k));
        n1 = dim(E1);
        E2 = ellipsoid.generateRandom(bools(2,k));
        n2 = dim(E2);
        E = cartProd(E1,E2);
        Y1 = randPoint(E1,n1);
        Y2 = randPoint(E2,n2);
        Y = combineVec(Y1,Y2);
        if ~in(E,Y)
            res = false;
            break;
        end
        % "in" check only works if at least one ellipsoid is non-degenerate
        if ((~E.isdegenerate || ~E1.isdegenerate) && ~in(project(E,1:n1),E1)) ||...
           ((~E.isdegenerate || ~E2.isdegenerate) && ~in(project(E,n1+1:n1+n2),E2))
            res = false;
            break;
        end
    end
    if ~res
        break;
    end
end
if ~res
    disp('testLongDuration_ellipsoid_cartProd failed');
else
    disp('testLongDuration_ellipsoid_cartProd successful');
end

%------------- END OF CODE --------------