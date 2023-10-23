function res = testLong_ellipsoid_radius
% testLong_ellipsoid_radius - unit test function of radius
%
% Syntax:
%    res = testLong_ellipsoid_radius
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       19-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% random tests
res = true;
nrOfTests = 100;
bools = [false,true];
for i=1:nrOfTests
    
    % random dimension
    n = randi(30);
    for k=1:2
        %%% generate all variables necessary to replicate results
        E = ellipsoid.generateRandom('Dimension',n,'IsDegenerate',bools(k));
        %%%
        E = ellipsoid(E.Q);
        [U,~,~] = svd(E.Q);
        IntE = interval(U'*E);
        R = rad(IntE);

        for j=1:dim(E)
            if ~all(withinTol(radius(E,j),R(1:j),E.TOL))
                res = false; return;
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
