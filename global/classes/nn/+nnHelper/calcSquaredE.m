function E_quad = calcSquaredE(E1, E2, isEqual)
% calcSquaredE - computes the multiplicative E1' * E2
%
% Syntax:
%    E_quad = nnHelper.calcSquaredE(E1, E2, isEqual)
%
% Inputs:
%    E1 - exponential matrix 1
%    E2 - exponential matrix 2
%    isEqual - whether they are equal for optimizations
%
% Outputs:
%    E_quad = result of E1' * E2 as row vector
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/quadMap, nnHelper/calcSquared

% Authors:       Tobias Ladner
% Written:       28-March-2022
% Last update:   23-June-2022 (performance optimizations)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% E1, E2 exponential matrices; calculate E1'*E2
if ~isempty(E1) && ~isempty(E2)
    [G1_ind, G2_ind, G1_ind2, G2_ind2] = nnHelper.calcSquaredGInd(E1(1, :), E2(1, :), isEqual);
    E_quad = [E1(:, G1_ind) + E2(:, G2_ind), E1(:, G1_ind2) + E2(:, G2_ind2)];

else
    E_quad = [];
end
end

% ------------------------------ END OF CODE ------------------------------
