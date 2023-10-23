function res = isemptyobject(hyp)
% isemptyobject - checks whether a constrained hyperplane contains any
%    information at all; consequently, the set is interpreted as the empty
%    set 
%
% Syntax:
%    res = isemptyobject(hyp)
%
% Inputs:
%    hyp - conHyperplane object
%
% Outputs:
%    res - true/false
%
% Example: 
%    hyp = conHyperplane([1;-1],1);
%    isemptyobject(hyp); % false
%    hyp = conHyperplane();
%    isemptyobject(hyp); % true
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

res = false(size(hyp));
% loop over class-arrays
for i=1:size(hyp,1)
    for j=1:size(hyp,2)
        res(i,j) = aux_checkIfEmpty(hyp(i,j));
    end
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkIfEmpty(hyp)

    res = isa(hyp.h,'halfspace') ...
        && isnumeric(hyp.h.c) && isempty(hyp.h.c) ...
        && isnumeric(hyp.h.d) && isscalar(hyp.h.d) && hyp.h.d == 0 ...
        && isnumeric(hyp.C) && isempty(hyp.C) ...
        && isnumeric(hyp.d) && isscalar(hyp.d) && hyp.d == 0;

end

% ------------------------------ END OF CODE ------------------------------
