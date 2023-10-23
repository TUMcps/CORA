function zB = quadMap(zB,Q)
% quadMap - computes \{Q_{ijk}*x_j*x_k|x \in Z\}
%
% Syntax:
%    zB = quadMap(zB,Q)
%
% Inputs:
%    zB - zonoBundle object
%    Q - quadratic coefficients as a cell of matrices
%
% Outputs:
%    zB - zonoBundle object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/quadMap

% Authors:       Niklas Kochdumper
% Written:       13-June-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{zB,'att','zonoBundle'}; ...
                {Q,'att','cell'}});

for i=1:zB.parallelSets
    zB.Z{i} = quadMap(zB.Z{i},Q);
end

% ------------------------------ END OF CODE ------------------------------
