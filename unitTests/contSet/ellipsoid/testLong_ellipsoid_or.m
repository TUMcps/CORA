function res = testLong_ellipsoid_or
% testLong_ellipsoid_or - unit test function of or
%
% Syntax:
%    res = testLong_ellipsoid_or
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
% Written:       14-October-2019
% Last update:   16-March-2021
% Last revision: 22-September-2024 (MW, integrate components)

% ------------------------------ BEGIN CODE -------------------------------

% double
assert(aux_orPoint());

% ellipsoid
assert(aux_orEllipsoid());

% test completed
res = true;

end


% Auxiliary functions -----------------------------------------------------

function res = aux_orPoint()

res = true;
nRuns = 2;
bools = [false,true];
% smaller dims since halfspaces and vertices are involved
for i=[2,5:5:10]
    for j=1:nRuns
        for k=1:2 
            try
                %%% generate all variables necessary to replicate results
                E = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',bools(k));
                V = randn(dim(E),2*dim(E));
                %%%
                
                for m=1:2
                    [U,S,V] = svd(V);
                    S(1,1) = bools(m)*S(1,1);
                    V = U*S*V';
                    % Eo = or(E,V);
                    Eo = E;
                    for ii=1:size(V,2)
                        Eo = or(Eo,V(:,ii));
                    end
                    
                    % check whether V are contained in Eo
                    assertLoop(all(contains(Eo,V)),i,j,k,m)
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

end

function res = aux_orEllipsoid()
% run test for ellipsoid/or

res = true;
nRuns = 2;
bools = [false,true];
for i=10:5:15
    for j=1:nRuns
        for k=1:2
            try
                E1 = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',false);
                E2 = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',false);
                E3 = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',bools(k)); 
                % compute outer approx (check overloaded syntax)
                Eo1 = E1 | E2;
                % check regular syntax
                Eo2 = or(E1,{E2,E3},'outer');
                % inner not implemented
%                 % compute inner approx
%                 Ei1 = or(E1,E2,'inner');
%                 Ei2 = or(E1,{E2,E3},'inner');
%                 % check if inner contained in outer
%                 if ~contains(Eo1,Ei1) || ~contains(Eo2,Ei2)
%                     res = false;
%                     break;
%                 end
                Y1 = randPoint(E1,i,'extreme');
                Y2 = randPoint(E2,i,'extreme');
                Y3 = randPoint(E3,i,'extreme');
                assertLoop(all(contains(Eo1,Y1)),i,j,k)
                assertLoop(all(contains(Eo1,Y2)),i,j,k)
                assertLoop(all(contains(Eo2,Y1)),i,j,k)
                assertLoop(all(contains(Eo2,Y2)),i,j,k)
                assertLoop(all(contains(Eo2,Y3)),i,j,k)

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

end

% ------------------------------ END OF CODE ------------------------------
