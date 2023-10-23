function res = testLong_component_ellipsoid_minkDiffEllipsoid
% testLong_component_ellipsoid_minkDiffEllipsoid - unit test
%
% Syntax:
%    res = testLong_component_ellipsoid_minkDiffEllipsoid
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
% Written:       17-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
nRuns = 2;
for i=10:5:15
    for j=1:nRuns
        try
            % for outer approx, both need to be degenerate, otherwise "and"
            % will throw error 
            %%% generate all variables necessary to replicate results
            E1 = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',false);
            E2 = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',false);
            m = ceil(dim(E2)*rand);
            Y1 = randPoint(E1,dim(E1),'extreme');
            Y2 = randPoint(E2,dim(E2),'extreme');
            %%%return
            % test "direct syntax
            E1o = minkDiff(E1,E2);
            E1i = minkDiff(E1,E2,'inner');
            [U,S,~] = svd(E2.Q);
            S(m,m) = 0;
            % generate degenerate ellipsoid contained in E2
            E3 = ellipsoid(U*S*U',E2.q);
            E2i = minkDiff(E1,E3,'inner');
            if representsa_(E1o,'emptySet',eps)
                if ~(E1==E2) && contains(E1,E2)
                    res = false;
                    return;
                end
                continue;
            end
            if ~contains(E1o,E1i) || ~contains(E1o,E2i)
                res = false;
                return;
            end
            E3i = minus(E2,E3,'inner');
            if representsa_(E3i,'emptySet',eps)
                res = false;
                return;
            end
            % test cell array syntax
            E2o = minkDiff(E1,E2);
            if ~(E1o==E2o)
                res = false;
                return;
            end
            % check points
            Y  = sumPoints(Y1,-Y2);
            %check if Y \in Eo
            if ~contains(E1o,Y)
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

% ------------------------------ END OF CODE ------------------------------
