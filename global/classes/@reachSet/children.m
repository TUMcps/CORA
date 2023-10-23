function out = children(R,parent)
% children - return a list of indices of the children of this parent node
%
% Syntax:
%    ch = children(R,5)
%
% Inputs:
%    R - reachSet object
%    parent - index of the parent node
%
% Outputs:
%    out - list of children node indices
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet

% Authors:       Benedikt Seidl
% Written:       19-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

out = [];

for i = 1:length(R)
    if R(i).parent == parent
        out(end+1) = i;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
