function res = testLong_component_ellipsoid_distanceEllipsoid
% testLong_component_ellipsoid_distanceEllipsoid - unit test 
%    function of distanceEllipsoid
%
% Syntax:
%    res = testLong_component_ellipsoid_distanceEllipsoid
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
% Written:       18-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
nRuns = 2;
bools = [false,true];
for i=5:5:10
    for j=1:nRuns
        for k=1:2 
            try
                %%% generate all variables necessary to replicate results
                E1 = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',false);
                E2 = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',bools(k));
                N = 2*i;
                m = ceil(N*rand);
                B = randPoint(E1,N,'extreme');
                fac_s = rand;
                r = ceil(i*rand);
                k = k;
                %%%
                b = B(:,m);
                s = E1.q + fac_s*(b-E1.q);
                % E1 and E2 are intersecting
                E2 = ellipsoid(E2.Q,s);
                if distance(E1,E2)>E1.TOL
                    res = false;
                    return;
                end
                IntE1 = interval(E1);
                q3 = E1.q+2.1*rad(IntE1);
                [U,S,~] = svd(E1.Q);
                S(r,r) = bools(k)*S(r,r);
                % guaranteed to not intersect with E1
                E3 = ellipsoid(U*S*U',q3);
                if distance(E1,E3)<=E1.TOL
                    res = false;
                    return;
                end
            catch ME
                if strcmp(ME.identifier,'CORA:solverIssue')
                    disp('Randomly generated ellipsoids caused solver issues! Ignoring...');
                    continue;
                end
                rethrow(ME);
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
