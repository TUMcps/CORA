function G_ = squareGenMat(G)
% squareGenMat - compute the generator matrix for the squared 
%    polynomial zonotope
%
% Syntax:
%    G = nnHelper.squareGenMat(G)
%
% Inputs:
%    G - generator matrix
%
% Outputs:
%    G - squared generator matrix
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

if ~isempty(G)
    n = length(G);
    G_ = zeros(1, 0.5*n*(n + 1));
    temp = G' * G;
    cnt = n;

    for i = 1:n - 1
        G_(i) = temp(i, i);
        G_(cnt+1:cnt+n-i) = 2 * temp(i, i+1:n);
        cnt = cnt + n - i;
    end
    G_(n) = temp(end, end);
else
    G_ = G;
end
end

%------------- END OF CODE --------------