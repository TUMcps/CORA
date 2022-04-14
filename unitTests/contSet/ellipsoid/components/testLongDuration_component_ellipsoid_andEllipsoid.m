function res = test_ellipsoid_andEllipsoid
% testLongDuration_ellipsoid_andEllipsoid - unit test function of testLongDuration_ellipsoid_andEllipsoid
%
% Syntax:  
%    res = testLongDuration_ellipsoid_andEllipsoid
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
                E1 = ellipsoid.generateRandom(false,i);
                N = 10*i;
                b = randPoint(ellipsoid(E1.Q),N,'extreme');
                % use boundary points as well as interior points
                s = rand(1,N).*b;
                Y = E1.q + [b,s];
                % use point from Y as center to make sure intersection is not
                % empty
                E2 = ellipsoid.generateRandom(bools(k),i);
                m = ceil(N*rand);
                E2 = ellipsoid(E2.Q,Y(:,m));
                Eo = E1&E2;
                if isempty(Eo)
                    res = false;
                    break;
                end

                % Check that any point of Y that is both in E1 and E2 is
                % also in the intersection, i.e., Eo
                for i_y = 1:size(Y,2)
                    if in(E1,Y(:,i_y)) && in(E2,Y(:,i_y)) && ~in(Eo,Y(:,i_y))
                        res = false;
                        break;
                    end
                end


                % only works if E1 and E2 full-dimensional, otherwise Eo and Ei
                % are both degenerate
                if ~bools(k)
                    Ei = and(E1,E2,'i');
                    if ~in(Eo,Ei)
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
%------------- END OF CODE --------------