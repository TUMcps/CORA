function res = isemptyobject(ls)
% isemptyobject - checks whether a level set contains any information at
%    all; consequently, the set is interpreted as the empty set 
%
% Syntax:
%    res = isemptyobject(ls)
%
% Inputs:
%    ls - levelSet object
%
% Outputs:
%    res - true/false
%
% Example: 
%    syms x y; eq = x^2 - y;
%    ls = levelSet(eq,[x;y],'<=');
%    isemptyobject(ls); % false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       24-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = false(size(ls));
% loop over class-arrays
for i=1:size(ls,1)
    for j=1:size(ls,2)
        res(i,j) = aux_checkIfEmpty(ls(i,j));
    end
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkIfEmpty(ls)

    res = isnumeric(ls.eq) && isempty(ls.eq) ...
        && isnumeric(ls.vars) && isempty(ls.vars) ...
        && isnumeric(ls.compOp) && isempty(ls.compOp) ...
        && isnumeric(ls.funHan) && isempty(ls.funHan) ...
        && isnumeric(ls.der) && isempty(ls.der) ...
        && isnumeric(ls.dim) && (isempty(ls.dim) || ls.dim == 0) ...
        && isnumeric(ls.solved) && isempty(ls.solved) ...
        && islogical(ls.solvable) && ~ls.solvable;

end

% ------------------------------ END OF CODE ------------------------------
