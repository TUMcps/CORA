function res = testLongDuration_component_ellipsoid_andEllipsoid
% testLongDuration_component_ellipsoid_andEllipsoid - unit test function of
%    testLongDuration_component_ellipsoid_andEllipsoid
%
% Syntax:  
%    res = testLongDuration_component_ellipsoid_andEllipsoid
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

% Author:       Victor Gassmann
% Written:      17-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
nRuns = 2;
bools = [false,true];
for i=10:5:15
    for j=1:nRuns
        for k=1:2
            try
                %%% generate all variables necessary to replicate results
                N = 10*i;
                E1 = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',false);
                E2 = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',bools(k));
                % to avoid false results due to numerical issues, do not
                % use "exact" boundary points
                b = 0.99*randPoint(ellipsoid(E1.Q),N,'extreme');
                fac_b = rand(1,N);
                m = ceil(N*rand);
                %%%
                % use boundary points as well as interior points
                s = fac_b.*b;
                Y = E1.q + [b,s];
                % use point from Y as center to make sure intersection is not
                % empty
                E2 = ellipsoid(E2.Q,Y(:,m));
                Eo = E1&E2;
                if isempty(Eo)
                    res = false;
                    break;
                end

                % Check that any point of Y that is both in E1 and E2 is
                % also in the intersection, i.e., Eo
                for i_y = 1:size(Y,2)
                    if contains(E1,Y(:,i_y)) && contains(E2,Y(:,i_y)) && ~contains(Eo,Y(:,i_y))
                        res = false;
                        break;
                    end
                end


                % only works if E1 and E2 full-dimensional, otherwise Eo and Ei
                % are both degenerate
                if ~bools(k)
                    Ei = and(E1,E2,'inner');
                    if ~contains(Eo,Ei)
                        res = false;
                        break;
                    end
                end
            catch ME
            if strcmp(ME.identifier,'CORA:solverIssue')
                disp('Randomly generated ellipsoids caused solver issues! Ignoring...');
                continue;
            end
            rethrow(ME);
            end
        end
        if ~res
            break;
        end
    end
    if ~res
        break;
    end
end

% save if failed
if ~res
    path = pathFailedTests(mfilename());
    save(path,'E1','E2','b','N','fac_b','m');
end
%------------- END OF CODE --------------