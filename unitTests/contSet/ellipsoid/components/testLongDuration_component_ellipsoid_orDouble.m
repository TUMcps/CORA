function res = testLongDuration_component_ellipsoid_orDouble
% testLongDuration_component_ellipsoid_orDouble - unit test function of 
% ellipsoid_orDouble
%
% Syntax:  
%    res = testLongDuration_component_ellipsoid_orDouble
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
% smaller dims since halfspaces and vertices are involved
for i=5:5:10
    for j=1:nRuns
        for k=1:2 
            try
                E = ellipsoid.generateRandom(i,bools(k));
                % generate points randomly
                N = 2*dim(E);
                V = randn(dim(E),N);
                for m=1:2
                    [U,S,V] = svd(V);
                    S(1,1) = bools(m)*S(1,1);
                    V = U*S*V';
                    V_c = mat2cell(V,E.dim,ones(1,N));
                    Eo = or(E,V_c);
                    Ei = or(E,V_c,'i');
                    % check whether subset property holds
                    if ~in(Eo,Ei)
                        res = false;
                        break;
                    end
                    % check whether V are contained in Eo
                    if ~in(Eo,V)
                        res = false;
                        break;
                    end
                end
                if ~res
                    break;
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