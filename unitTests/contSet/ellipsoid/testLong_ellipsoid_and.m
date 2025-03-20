function res = testLong_ellipsoid_and
% testLong_ellipsoid_and - unit test function of intersection
%
% Syntax:
%    res = testLong_ellipsoid_and
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
% Written:       17-October-2019
% Last update:   17-March-2021
% Last revision: 22-September-2024 (MW, integrate component tests)

% ------------------------------ BEGIN CODE -------------------------------

% intersection with ellipsoids
assert(aux_andEllipsoid());

% intersection with polytope
assert(aux_andPolytope());

% intersection with hyperplanes
assert(aux_andHyperplane());

% combine results
res = true;

end


% Auxiliary functions -----------------------------------------------------

function res = aux_andEllipsoid()

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
                assertLoop(~representsa_(Eo,'emptySet',eps),i,j,k)

                % Check that any point of Y that is both in E1 and E2 is
                % also in the intersection, i.e., Eo
                for i_y = 1:size(Y,2)
                    if contains(E1,Y(:,i_y)) && contains(E2,Y(:,i_y))
                        assertLoop(contains(Eo,Y(:,i_y)),i,j,k)
                    end
                end


                % only works if E1 and E2 full-dimensional, otherwise Eo and Ei
                % are both degenerate
                if ~bools(k)
                    Ei = and_(E1,E2,'inner');
                    assertLoop(contains(Eo,Ei),i,j,k)
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

function res = aux_andPolytope()
% compute polytope and

res = true;
nRuns = 2;
bools = [false,true];
for i=10:5:15
    for j=1:nRuns
        for k=1:2
            %%% generate all variables necessary to replicate results
            E = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',bools(k));
            c = randn(E.dim,1);
            B = randPoint(E,2*dim(E),'extreme');
            m = ceil(2*dim(E)*rand);
            fac_s = rand;
            v = randn(dim(E),1);
            %%%
            % normalize tangent
            c = c/norm(c);
            [val,x] = supportFunc(E,c);
            % completely contains E
            P1 = polytope(c',val);
            assertLoop((E & P1) == E,i,j)

            % completely outside (ignoring touching)
            P2 = polytope(-c',-val);
            res2 = E & P2;
   
            assertLoop(1/cond(res2.Q) <= E.TOL,i,j,k)
            assertLoop(all(withinTol(res2.q,x,E.TOL)),i,j,k)

            % choose random point within E
            s = E.q + fac_s*(B(:,m)-E.q);
            P3 = polytope(v',v'*s);
            assertLoop(~representsa_(E & P3,'emptySet',eps),i,j,k)
        end
    end
end

end

function res = aux_andHyperplane()
% compute hyperplane and

res = true;
nRuns = 2;
bools = [false,true];
for i=10:5:15
    for j=1:nRuns
        for k=1:2
            E = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',bools(k));
            % construct hyperplane that intersects and check
            % select random boundary point
            N = 2*E.dim;
            m = ceil(N*rand);
            B = randPoint(E,N,'extreme');
            b = B(:,m);
            % select sample somewhere between b and E.q
            s = E.q + rand*(b-E.q);
            % generate random direction
            d = randn(E.dim,1);
            d = d/norm(d);
            % generate hyperplane
            H = polytope([],[],d',d'*s);
            assertLoop(~representsa_(E&H,'emptySet',eps),i,j,k)
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
