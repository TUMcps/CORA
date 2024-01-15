function res = isemptyobject(hs)
% isemptyobject - checks whether a halfspace contains any information at
%    all; consequently, the set is interpreted as the empty set 
%
% Syntax:
%    res = isemptyobject(hs)
%
% Inputs:
%    hs - halfspace object
%
% Outputs:
%    res - true/false
%
% Example: 
%    hs = halfspace([1;-1],1);
%    isemptyobject(hs); % false
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

res = false(size(hs));
% loop over class-arrays
for i=1:size(hs,1)
    for j=1:size(hs,2)
        res(i,j) = aux_checkIfEmpty(hs(i,j));
    end
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkIfEmpty(hs)

    res = isnumeric(hs.c) && isempty(hs.c) ...
        && isnumeric(hs.d) && isscalar(hs.d) && hs.d == 0;

end

% ------------------------------ END OF CODE ------------------------------
