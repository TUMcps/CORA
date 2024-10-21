function res = testLong_ellipsoid_minkDiff
% testLong_ellipsoid_minkDiff - unit test function of minkDiff
%
% Syntax:
%    res = testLong_ellipsoid_minkDiff
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
% See also: none

% Authors:       Victor Gassmann
% Written:       15-March-2021
% Last update:   ---
% Last revision: 22-September-2024 (MW, integrate components)

% ------------------------------ BEGIN CODE -------------------------------

% point
assert(aux_minkDiffPoint());

% ellipsoid
assert(aux_minkDiffEllipsoid());


% combine results
res = true;

end


% Auxiliary functions -----------------------------------------------------

function res = aux_minkDiffPoint()

res = true;
nRuns = 5;
bools = [false,true];

% smaller dims since halfspaces and vertices are involved
for i=10:5:15
    for j=1:nRuns

        % degenerate and non-degenerate
        for k=1:2 
            % init random ellipsoid
            E = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',bools(k));
            % generate points randomly
            V = randn(dim(E),2*dim(E));
            
            for m=1:2
                [U,S,V] = svd(V);
                S(1,1) = bools(m)*S(1,1);
                V = U*S*V';
                Eo = E; Ei = E;
                for ii=1:size(V,2)
                    Eo = minkDiff(Eo,V(:,ii));
                    Ei = minkDiff(Ei,V(:,ii),'inner');
                end
                Eres = ellipsoid(E.Q,E.q-sum(V,2));
                % check if the same
                assertLoop((Eres==Eo),i,j,k,m)
                assertLoop((Eres==Ei),i,j,k,m)
            end
        end

    end
end

end

function res = aux_minkDiffEllipsoid()

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
                assertLoop((E1==E2) || ~contains(E1,E2),i,j)
                continue;
            end
            assertLoop(contains(E1o,E1i),i,j)
            assertLoop(contains(E1o,E2i),i,j)

            E3i = minus(E2,E3,'inner');
            assertLoop(~representsa_(E3i,'emptySet',eps),i,j)

            % test cell array syntax
            E2o = minkDiff(E1,E2);
            assertLoop((E1o==E2o),i,j)

            % check points
            Y  = sumPoints(Y1,-Y2);
            %check if Y \in Eo
            assertLoop(contains(E1o,Y),i,j)

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
