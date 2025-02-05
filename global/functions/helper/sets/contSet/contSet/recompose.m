function S_out = recompose(S)
% recompose - recompose a cell array of sets via the Cartesian product to
%    one set of higher dimension
%
% Syntax:
%    S_out = recompose(S)
%
% Inputs:
%    S - cell array of contSet objects/numeric, numeric
%
% Outputs:
%    S - contSet object
%
% Example:
%    S = {zonotope(0),zonotope([1;-1])};
%    S_out = recompose(S);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/recompose, decompose, linearSys/priv_reach_decomp

% Authors:       Mark Wetzlinger
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% single vector
if isnumeric(S)
    S_out = S;
    return
end

% concatenate vector
if all(cellfun(@isnumeric,S,'UniformOutput',true))
    S_out = vertcat(S{:});
    return
end

% sequentially compute the Cartesian product
% for safety reasons, we read out the first element (unknown whether
% all contSet/cartProd_ functions support 'empty numeric' case)
S_out = S{1};
for i=2:numel(S)
    S_out = cartProd(S_out,S{i});
end

% ------------------------------ END OF CODE ------------------------------
