function res = testLong_linearSys_falsify
% testLong_linearSys_falsify - unit test for falsification
%
% Syntax:
%    res = testLong_linearSys_falsify
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

% Authors:       Niklas Kochdumper
% Written:       27-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    alg = {'singleOpt','multiOpt','monteCarlo'};

    for n = [1 2 3]
        for m = [0 1 2]
            for p = [1 2]

                % create random linear system
                B = []; D = [];
                if m > 0
                    B = rand(n,m); D = rand(p,m);
                end

                sys = linearSys(rand(n),B,rand(n,1),rand(p,n),D,rand(p,1));

                % create random reachability problem
                params = [];
                params.R0 = aux_randomSet(n);

                if m > 0
                    params.U = aux_randomSet(m);
                end

                params.tFinal = rand();

                % create random specification
                P = aux_randomSet(p);

                if randi(2) == 1
                    spec = specification(P,'unsafeSet');
                else
                    spec = specification(P,'safeSet');
                end

                % loop over all falsification algorithms
                for k = 1:length(alg)

                    options.falsifyAlg = alg{k};
                    options.maxTime = 2;

                    falsify(sys,params,options,spec);
                end
            end
        end
    end

    % test successfull if it runs through without errors
    res = true;

end


% Auxiliary functions -----------------------------------------------------

function S = aux_randomSet(dims)
% generate a random set, where the set representation is also picked at
% random

    type = randi(5);

    switch type
        case 1
            S = interval.generateRandom('Dimension',dims);
        case 2
            S = zonotope.generateRandom('Dimension',dims);
        case 3
            S = polytope.generateRandom('Dimension',dims);
        case 4
            S = conZonotope.generateRandom('Dimension',dims);
        case 5
            S = zonoBundle.generateRandom('Dimension',dims);
    end
end

% ------------------------------ END OF CODE ------------------------------
