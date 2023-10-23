function [G_quad] = calcSquaredG(G1, G2, isEqual)
% calcSquaredG - computes the multiplicative G1' * G2
%
% Syntax:
%    G_quad = nnHelper.calcSquaredG(G1, G2, isEqual)
%    G_quad = nnHelper.calcSquaredG(G1, G2)
%
% Inputs:
%    G1 - row vector of generator matrix 1
%    G2 - row vector of generator matrix 2
%    isEqual - whether they are equal for optimizations
%
% Outputs:
%    G_quad = result of G1' * G2 as row vector
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/quadMap, nnHelper/calcSquared

% Authors:       Tobias Ladner
% Written:       28-March-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 3
    isEqual = false;
end

if ~isempty(G1) && ~isempty(G2)
    if isEqual
        temp = G1' * G2;

        % we can ignore the left lower triangle in this case
        % as it's the same as the right upper triangle
        % -> double right upper triangle
        n = length(G1);
        G_quad = zeros(1, 0.5*n*(n + 1));
        cnt = n;

        for i = 1:n - 1
            G_quad(i) = temp(i, i);
            G_quad(cnt+1:cnt+n-i) = 2 * temp(i, i+1:n);
            cnt = cnt + n - i;
        end
        G_quad(n) = temp(end, end);
    else
        % calculate all values
        G_quad = G1' * G2;
        G_quad = reshape(G_quad, 1, []); % row vector
    end
else
    G_quad = [];
end
end

% ------------------------------ END OF CODE ------------------------------
