function res = isemptyobject(simRes)
% isemptyobject - checks if a simResult object is empty
%
% Syntax:
%    res = isemptyobject(simRes)
%
% Inputs:
%    simRes - simResult object
%
% Outputs:
%    res - true/false
%
% Example:
%    simRes = simResult();
%    isemptyobject(simRes)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       01-May-2023
% Last update:   22-May-2023 (MW, extend to class arrays)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[r,c] = size(simRes);
res = false(r,c);

% loop over all objects
for i=1:r
    for j=1:c
        % check time
        res(i,j) = isempty(simRes(i,j).t);
    end
end

% ------------------------------ END OF CODE ------------------------------
