function res = isemptyobject(pHA)
% isemptyobject - checks if a parallelHybridAutomaton object is empty
%
% Syntax:
%    res = isemptyobject(pHA)
%
% Inputs:
%    pHA - parallelHybridAutomaton object
%
% Outputs:
%    res - true/false
%
% Example: 
%    res = isemptyobject(parallelHybridAutomaton());
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       19-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[r,c] = size(pHA);
res = false(r,c);

% loop over all locations
for i=1:r
    for j=1:c
        res(i,j) = isempty(pHA(i,j).bindsInputs);
    end
end

% ------------------------------ END OF CODE ------------------------------
