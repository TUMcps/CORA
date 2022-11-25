function res = testLongDuration_ellipsoid_minusDouble
% testLongDuration_ellipsoid_minusDouble - unit test function of testLongDuration_ellipsoid_minusDouble
%
% Syntax:  
%    res = testLongDuration_ellipsoid_minusDouble
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
% Written:      18-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
nRuns = 5;
bools = [false,true];
% smaller dims since halfspaces and vertices are involved
for i=10:5:15
    for j=1:nRuns
        for k=1:2 
            E = ellipsoid.generateRandom(i,bools(k));
            % generate points randomly
            N = 2*dim(E);
            V = randn(dim(E),N);
            for m=1:2
                [U,S,V] = svd(V);
                S(1,1) = bools(m)*S(1,1);
                V = U*S*V';
                V_c = mat2cell(V,E.dim,ones(1,N));
                Eo = minus(E,V_c);
                Ei = minus(E,V_c,'i');
                Eres = ellipsoid(E.Q,E.q-sum(V,2));
                % check if the same
                if ~(Eres==Eo) || ~(Eres==Ei)
                    res = false;
                    break;
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
    if ~res
        break;
    end
end
%------------- END OF CODE --------------