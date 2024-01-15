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
%    hyp = conHyperplane([1 -1],1);
%    isemptyobject(hyp); % false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       24-July-2023
% Last update:   09-January-2024 (MW, update to new constructor)
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

    res = isnumeric(hyp.a) && isempty(hyp.a) ...
        && isnumeric(hyp.b) && isempty(hyp.b) ...
        && isnumeric(hyp.C) && isempty(hyp.C) ...
        && isnumeric(hyp.d) && isempty(hyp.d);

end

% ------------------------------ END OF CODE ------------------------------
