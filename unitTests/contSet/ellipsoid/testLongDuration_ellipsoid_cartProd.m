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
        %%% generate all variables necessary to replicate results
        E1 = ellipsoid.generateRandom('IsDegenerate',bools(1,k));
        E2 = ellipsoid.generateRandom('IsDegenerate',bools(2,k));
        Y1 = randPoint(E1,dim(E1));
        Y2 = randPoint(E2,dim(E2));
        %%%
        n1 = dim(E1);
        n2 = dim(E2);
        E = cartProd(E1,E2);
        Y = combineVec(Y1,Y2);
        if ~contains(E,Y)
            res = false;
            break;
        end
        % "in" check only works if at least one ellipsoid is non-degenerate
        if ((isFullDim(E) || isFullDim(E1)) && ~contains(project(E,1:n1),E1)) ||...
           ((isFullDim(E) || isFullDim(E2)) && ~contains(project(E,n1+1:n1+n2),E2))
            res = false;
            break;
        end
    end
    if ~res
        break;
    end
end

if ~res
    path = pathFailedTests(mfilename());
    save(path,'E1','E2','Y1','Y2');
end

%------------- END OF CODE --------------