function [G1_ind, G2_ind, G1_ind2, G2_ind2] = calcSquaredGInd(G1, G2, isEqual)
% calcSquaredGInd - computes the indices of multiplicative G1' * G2
%
% Syntax:
%    [G1_ind, G2_ind, G1_ind2, G2_ind2] = nnHelper.calcSquaredGInd(G1, G2, isEqual)
%    G_quad = nnHelper.calcSquaredGInd(G1, G2)
%
% Inputs:
%    G1 - row vector of generator matrix 1
%    G2 - row vector of generator matrix 2
%    isEqual - whether they are equal for optimizations
%
% Outputs:
%    G1_ind = indices of G1 in the final matrix
%    G2_ind = indices of G2 in the final matrix
%    G1_ind2 = indices of G1 in the final matrix which occur twice
%    G2_ind2 = indices of G2 in the final matrix which occur twice
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nnHelper/calcSquaredG, nnHelper/calcSquaredE

% Authors:       Tobias Ladner
% Written:       23-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 3
    isEqual = false;
end

G1_ind = [];
G2_ind = [];
G1_ind2 = [];
G2_ind2 = [];

if ~isempty(G1) && ~isempty(G2)
    tempG1 = (1:length(G1))'.^(ones(size(G2)));
    tempG2 = (1:length(G2)).^(ones(size(G1))');

    if isEqual
        % we can ignore the left lower triangle in this case
        % as it's the same as the right upper triangle
        % -> double right upper triangle

        G1_ind = 1:length(G1);
        G2_ind = 1:length(G2);

        G1_ind2 = reshape(triu(tempG1, 1)', 1, []);
        G2_ind2 = reshape(triu(tempG2, 1)', 1, []);

        G1_ind2 = G1_ind2(G1_ind2 > 0);
        G2_ind2 = G2_ind2(G2_ind2 > 0);

    else
        % calculate all values
        G1_ind = tempG1;
        G2_ind = tempG2;

        G1_ind = reshape(G1_ind, 1, []);
        G2_ind = reshape(G2_ind, 1, []);
    end
end
end

% ------------------------------ END OF CODE ------------------------------
