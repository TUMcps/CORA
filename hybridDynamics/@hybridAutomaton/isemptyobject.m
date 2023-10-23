function res = isemptyobject(HA)
% isemptyobject - checks if a hybrid automaton object is empty
%
% Syntax:
%    res = isemptyobject(HA)
%
% Inputs:
%    HA - hybridAutomaton object
%
% Outputs:
%    res - true/false
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       16-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[r,c] = size(HA);
res = false(r,c);

% loop over all locations
for i=1:r
    for j=1:c
        res(i,j) = all(isemptyobject(HA(i,j).location));
    end
end

% ------------------------------ END OF CODE ------------------------------
