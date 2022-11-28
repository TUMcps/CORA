function G = cubeGenMat(G)
% cubeGenMat - compute the generator matrix for the cubed 
%    polynomial zonotope
%
% Syntax:
%    G = nnHelper.cubeGenMat(G)
%
% Inputs:
%    G - generator matrix
%
% Outputs:
%    G - cubed generator matrix
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:        Niklas Kochdumper, Tobias Ladner
% Written:       17-September-2021
% Last update:   ---
% Last revision: 28-March-2022 (TL)

%------------- BEGIN CODE --------------

G1 = [];
G2 = [];
G3 = [];

for i = 1:size(G, 2) - 1
    for j = i + 1:size(G, 2)
        G1 = [G1, 3 * G(i)^2 * G(j)];
    end
end

for i = 1:size(G, 2)
    for j = 1:i - 1
        G2 = [G2, 3 * G(i)^2 * G(j)];
    end
end

for i = 1:size(G, 2) - 1
    for j = i + 1:size(G, 2) - 1
        for k = j + 1:size(G, 2)
            G3 = [G3, 6 * G(i) * G(j) * G(k)];
        end
    end
end

G = [G.^3, G1, G2, G3];
end

%------------- END OF CODE --------------