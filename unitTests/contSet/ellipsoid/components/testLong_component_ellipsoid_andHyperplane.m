function res = testLong_component_ellipsoid_andHyperplane
% testLong_component_ellipsoid_andHyperplane - unit test function
%    of andHyperplane
%
% Syntax:
%    res = testLong_component_ellipsoid_andHyperplane
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
            H = conHyperplane(d',d'*s);
            if representsa_(E&H,'emptySet',eps)
                res = false;
                return;
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
